"""Unit tests for the MAICoS wrappers in fairmd.lipids.analib.maicos."""

import io
import json
import logging
import os
import shutil
import tempfile
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pytest
from numpy.testing import assert_allclose

pytestmark = [pytest.mark.sim1, pytest.mark.min]

TOY_SYSTEM = Path(__file__).parents[1] / "ToyData" / "Simulations.1" / "aa0"
TPR = TOY_SYSTEM / "md.tpr"
XTC = TOY_SYSTEM / "short.xtc"
LAST_ATOM = "C50"  # terminal methyl carbon of the Berger DPPC sn-1 tail
G3_ATOM = "C12"
# The trajectory has 121 frames at 1 ns spacing; skipping 118 ns keeps the last 3.
EQ_TIME = 118_000
N_FRAMES_AFTER_EQ = 3

requires_gmx = pytest.mark.skipif(shutil.which("gmx") is None, reason="GROMACS not installed")


@pytest.fixture
def universe():
    return mda.Universe(TPR, XTC)


@pytest.mark.parametrize(
    "obj",
    [np.arange(3.0), np.arange(6.0).reshape(3, 2), {"profile": np.arange(3.0)}],
    ids=["1d", "2d", "nested"],
)
def test_numpy_array_encoder_roundtrip(obj):
    from fairmd.lipids.analib.maicos import NumpyArrayEncoder

    buffer = io.StringIO()
    json.dump(obj, buffer, cls=NumpyArrayEncoder)
    decoded = json.loads(buffer.getvalue())

    expected = {k: v.tolist() for k, v in obj.items()} if isinstance(obj, dict) else obj.tolist()
    assert decoded == expected


@pytest.mark.parametrize(
    ("analysis", "kwargs", "columns"),
    [
        ("DensityPlanar", {"output": "out.json"}, {"out.json": ("profile", "dprofile", 1e3)}),
        ("DiporderPlanar", {"output": "out.json"}, {"out.json": ("profile", "dprofile", 1)}),
        (
            "DielectricPlanar",
            {"output_prefix": "out", "temperature": 323},
            {"out_perp.json": ("eps_perp", "deps_perp", 1), "out_par.json": ("eps_par", "deps_par", 1)},
        ),
    ],
)
@pytest.mark.filterwarnings("ignore:Your trajectory is too short:UserWarning")
def test_save_writes_profiles_in_nm(universe, tmp_path, monkeypatch, analysis, kwargs, columns):
    """Smoke-run each wrapper and check save() writes [z/nm, value, error] rows."""
    import fairmd.lipids.analib.maicos

    monkeypatch.chdir(tmp_path)
    water = universe.select_atoms("resname SOL")
    result = getattr(fairmd.lipids.analib.maicos, analysis)(water, **kwargs).run(stop=3)
    result.save()

    for filename, (value, error, scale) in columns.items():
        saved = json.loads((tmp_path / filename).read_text())
        expected = np.column_stack(
            [result.results.bin_pos / 10, scale * result.results[value], scale * result.results[error]]
        )
        # CompactJSONEncoder writes floats with 6 significant digits
        assert_allclose(saved, expected, rtol=1e-5)


@pytest.mark.parametrize(
    ("system", "expected"),
    [
        ({"TYPEOFSYSTEM": "miscellaneous"}, False),
        ({"TYPEOFSYSTEM": "lipid bilayer"}, True),
        ({"TYPEOFSYSTEM": "lipid bilayer", "WARNINGS": {"ORIENTATION": "warning"}}, False),
        ({"TYPEOFSYSTEM": "lipid bilayer", "WARNINGS": {"PBC": "hexagonal-box"}}, False),
        ({"TYPEOFSYSTEM": "lipid bilayer", "WARNINGS": {"NOWATER": True}}, False),
        ({"TYPEOFSYSTEM": "lipid bilayer", "WARNINGS": {"PBC": "rectangular-box"}}, True),
        ({"TYPEOFSYSTEM": "lipid bilayer", "WARNINGS": {}}, True),
        ({"TYPEOFSYSTEM": "lipid bilayer", "WARNINGS": None}, True),
    ],
)
def test_maicos_suitable(system, expected):
    from fairmd.lipids.analib.maicos import is_system_suitable_4_maicos

    assert is_system_suitable_4_maicos(system) is expected


def test_maicos_what_to_compute(caplog, logger):
    from fairmd.lipids import FMDL_SIMU_PATH, RCODE_ERROR
    from fairmd.lipids.analyze import computeMAICOS
    from fairmd.lipids.core import System

    with tempfile.TemporaryDirectory(dir=FMDL_SIMU_PATH) as path:
        for filename in ("file.xtc", "file.tpr", "file.pdb"):
            open(os.path.join(path, filename), "w").close()

        system = System(
            data={
                "DOI": "00.0000/abcd",
                "path": os.path.basename(path),
                "TYPEOFSYSTEM": "lipid bilayer",
                "SOFTWARE": "GROMACS",
                "TPR": [["file.tpr"]],
                "TRJ": [["file.xtc"]],
                "ID": 9999,
                "COMPOSITION": {
                    "SOL": {"NAME": "SPC", "MAPPING": "mappingSPCwater.yaml", "COUNT": 1000},
                    "DPPC": {"NAME": "DPPC", "MAPPING": "mappingDPPCberger.yaml", "COUNT": 100},
                },
            },
        )
        with caplog.at_level(logging.INFO):
            rcode = computeMAICOS(system, logging.getLogger("test_logger"), ffonly=False)
        line = next(line for line in caplog.text.splitlines() if "Files to be computed:" in line)
        assert "TotalDensity" in line
        assert "DiporderWater" in line
        assert rcode == RCODE_ERROR


def _center_mda(universe, path, **kwargs):
    from fairmd.lipids.analib.maicos import traj_centering_for_maicos_mda

    return traj_centering_for_maicos_mda(universe, path, LAST_ATOM, EQ_TIME, **kwargs)


def _center_mda_parallel(universe, path, **kwargs):
    from fairmd.lipids.analib.maicos import traj_centering_for_maicos_mda_parallel

    return traj_centering_for_maicos_mda_parallel(universe, path, LAST_ATOM, EQ_TIME, n_jobs=2, **kwargs)


def _center_gromacs(_universe, path, **kwargs):
    from fairmd.lipids.analib.maicos import traj_centering_for_maicos_gromacs

    return traj_centering_for_maicos_gromacs(path, str(XTC), str(TPR), LAST_ATOM, G3_ATOM, EQ_TIME, **kwargs)


@pytest.mark.parametrize(
    ("center", "output_name"),
    [
        (_center_mda, "whole.xtc"),
        (_center_mda_parallel, "whole.xtc"),
        (_center_gromacs, "centered.xtc"),
    ],
    ids=["mda", "mda_parallel", "gromacs"],
)
def test_centering_reuses_existing_output(universe, tmp_path, center, output_name):
    existing = tmp_path / output_name
    existing.write_text("dummy")

    assert center(universe, str(tmp_path)) == str(existing)
    assert existing.read_text() == "dummy"


def _bounding_box_center(atomgroup):
    return (atomgroup.positions.min(axis=0) + atomgroup.positions.max(axis=0)) / 2


@pytest.mark.parametrize(
    ("center", "reference_atom", "group_center"),
    [
        (_center_mda, LAST_ATOM, lambda atomgroup: atomgroup.center_of_mass()),
        # `gmx trjconv -center` centres the midpoint of the group's bounding box
        pytest.param(_center_gromacs, G3_ATOM, _bounding_box_center, marks=requires_gmx),
    ],
    ids=["mda", "gromacs"],
)
def test_centering_puts_reference_at_box_center(universe, tmp_path, center, reference_atom, group_center):
    centered = mda.Universe(TPR, center(universe, str(tmp_path), recompute=True))

    assert centered.trajectory.n_frames == N_FRAMES_AFTER_EQ
    reference = centered.select_atoms(f"name {reference_atom}")
    for ts in centered.trajectory:
        # 0.1 Å is well above the 0.01 Å XTC precision but far below any bilayer offset
        assert_allclose(group_center(reference)[2], ts.dimensions[2] / 2, atol=0.1)


def test_parallel_centering_matches_sequential(universe, tmp_path):
    sequential_dir = tmp_path / "sequential"
    parallel_dir = tmp_path / "parallel"
    sequential_dir.mkdir()
    parallel_dir.mkdir()

    sequential = mda.Universe(TPR, _center_mda(universe, str(sequential_dir)))
    parallel = mda.Universe(TPR, _center_mda_parallel(universe, str(parallel_dir)))

    assert parallel.trajectory.n_frames == sequential.trajectory.n_frames
    for _ in zip(sequential.trajectory, parallel.trajectory, strict=True):
        # XTC stores positions to 0.01 Å
        assert_allclose(parallel.atoms.positions, sequential.atoms.positions, atol=0.01)
