"""Unit tests for MAICoS-specific functionality."""

import json
import logging
import os
import shutil
import tempfile
import warnings

import pytest
import pytest_check as check

# run only on sim2 mocking data
pytestmark = [pytest.mark.sim1, pytest.mark.min]


@pytest.fixture
def tip4p_trajectory(tmpdir):
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.coordinates.memory import MemoryReader

    folder = str(tmpdir)
    pdb_content = """TITLE     Single 4-site water
CRYST1   63.646   63.646   86.762  90.00  90.00  90.00 P 1           1
ATOM      1  OW  wate    1       4.370  31.010  61.030  1.00  0.00
ATOM      2  HW1 wate    1       3.800  30.820  61.780  1.00  0.00
ATOM      3  HW2 wate    1       5.250  30.760  61.320  1.00  0.00
ATOM      4  MW  wate    1       4.410  30.950  61.170  1.00  0.00
END
"""
    tmpfile = os.path.join(folder, "test_water_maicos.pdb")
    with open(tmpfile, "w") as f:
        f.write(pdb_content)
    with pytest.warns(UserWarning, match="Element information is missing"):
        u = mda.Universe(tmpfile, format="PDB")
    nul_dim = u.dimensions
    trajectory = np.zeros((10, u.atoms.n_atoms, 3), dtype=np.float32)
    shift_v = np.array([0.0, 0.0, -0.5], dtype=np.float32)
    for i in range(10):
        trajectory[i] = u.atoms.positions.copy() + shift_v * i
    u.load_new(trajectory, order="fac", format=MemoryReader)
    for ts in u.trajectory:
        ts.dimensions = nul_dim
    yield u, folder
    shutil.rmtree(folder)


@pytest.fixture
def tip4p_random_rotation(tmpdir):
    """TIP4P water rotating continuously for 1000 1-ps frames in memory."""
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.coordinates.memory import MemoryReader

    folder = str(tmpdir)
    pdb_content = """TITLE     Single 4-site water
CRYST1   63.646   63.646   86.762  90.00  90.00  90.00 P 1           1
ATOM      1  OW  wate    1       4.370  31.010  61.030  1.00  0.00
ATOM      2  HW1 wate    1       3.800  30.820  61.780  1.00  0.00
ATOM      3  HW2 wate    1       5.250  30.760  61.320  1.00  0.00
ATOM      4  MW  wate    1       4.410  30.950  61.170  1.00  0.00
END
"""
    tmpfile = os.path.join(folder, "test_water_maicos.pdb")
    with open(tmpfile, "w") as f:
        f.write(pdb_content)
    with pytest.warns(UserWarning, match="Element information is missing"):
        u = mda.Universe(tmpfile, format="PDB")

    box = u.dimensions.copy()
    origin = u.atoms.positions[0].copy()
    relative = u.atoms.positions - origin
    rng = np.random.default_rng(seed=0)
    n_frames = 400
    trajectory = np.empty((n_frames, u.atoms.n_atoms, 3), dtype=np.float32)
    rotation = np.eye(3)
    for frame in range(n_frames):
        axis = rng.normal(size=3)
        axis /= np.linalg.norm(axis)
        angle = rng.uniform(-np.pi / 12, np.pi / 12)
        x, y, z = axis
        c, s = np.cos(angle), np.sin(angle)
        increment = np.array(
            [
                [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
                [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
                [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
            ],
        )
        rotation = increment @ rotation
        trajectory[frame] = origin + relative @ rotation.T

    u.load_new(trajectory, order="fac", format=MemoryReader, dt=1.0)
    for ts in u.trajectory:
        ts.dimensions = box
    yield u, folder
    shutil.rmtree(folder)


def test_density_planar(tip4p_trajectory, tmp_path):
    from fairmd.lipids.analib.maicos import DensityPlanar

    u, _ = tip4p_trajectory
    u.add_TopologyAttr("elements")
    u.atoms.elements = ["O", "H", "H", "Dummy"]
    ofname = tmp_path / "output.json"
    dp = DensityPlanar(
        u.atoms,
        dens="electron",
        zmin=0,
        zmax=20,
        unwrap=False,
        bin_width=1,
        pack=False,
        output=ofname,
    )
    dp.run()
    dp.save()
    assert os.path.isfile(ofname)

    with open(ofname) as output_file:
        output = json.load(output_file)
    assert len(output) == 20  # noqa: PLR2004
    assert all(len(row) == 3 for row in output)  # noqa: PLR2004
    assert max(row[1] for row in output) > 0
    check.equal(output[0][1], pytest.approx(0.0, abs=1e-6))
    check.equal(output[14][1], pytest.approx(0.493686, abs=1e-6))
    check.equal(output[14][2], pytest.approx(0.294148, abs=1e-6))
    check.equal(output[-1][1], pytest.approx(0.0, abs=1e-6))


def test_dielectric_planar(tip4p_random_rotation, tmp_path):
    from fairmd.lipids.analib.maicos import DielectricPlanar

    u, _ = tip4p_random_rotation
    u.add_TopologyAttr("resids", values=[1])
    u.add_TopologyAttr("resnames", values=["WAT"])
    u.add_TopologyAttr("charges")
    u.atoms.charges = [-0.848, 0.424, 0.424, 0.0]
    output_prefix = str(tmp_path / "dielectric")

    class IgnoreDielectricWarning(logging.Filter):
        def filter(self, record):
            return "Setting `zmin` and `zmax` might cut off molecules." not in record.getMessage()

    root_logger = logging.getLogger()
    warning_filter = IgnoreDielectricWarning()
    root_logger.addFilter(warning_filter)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"Your trajectory does not provide sufficient statistics to estimate a correlation time\.",
            category=UserWarning,
            module=r"maicos\.core\.base",
        )
        dielectric = DielectricPlanar(
            u.atoms,
            temperature=300,
            zmin=0,
            zmax=20,
            unwrap=False,
            bin_width=1,
            pack=False,
            output_prefix=output_prefix,
        )
        dielectric.run()
    root_logger.removeFilter(warning_filter)
    dielectric.save()

    expected_val_17 = {
        "_perp.json": (-0.0662769, 0.0137618),
        "_par.json": (0.0840072, 0.00411034),
    }
    for suffix in ("_perp.json", "_par.json"):
        output = f"{output_prefix}{suffix}"
        assert os.path.isfile(output)
        with open(output) as output_file:
            result = json.load(output_file)
        assert len(result) == 20  # noqa: PLR2004
        assert all(len(row) == 3 for row in result)  # noqa: PLR2004
        check.equal(result[0][0], pytest.approx(0.05, abs=1e-6))
        check.equal(result[17][1], pytest.approx(expected_val_17[suffix][0], rel=1e-4))
        check.equal(result[17][2], pytest.approx(expected_val_17[suffix][1], rel=1e-4))


def test_diporder_planar(tip4p_trajectory, tmp_path):
    from fairmd.lipids.analib.maicos import DiporderPlanar

    u, _ = tip4p_trajectory
    u.add_TopologyAttr("resids", values=[1])
    u.add_TopologyAttr("resnames", values=["WAT"])
    u.add_TopologyAttr("charges")
    u.atoms.charges = [-0.848, 0.424, 0.424, 0.0]
    output = tmp_path / "diporder.json"
    dp = DiporderPlanar(
        u.atoms,
        zmin=0,
        zmax=20,
        unwrap=False,
        bin_width=1,
        pack=False,
        output=output,
    )
    dp.run()
    dp.save()

    assert os.path.isfile(output)
    with open(output) as output_file:
        result = json.load(output_file)
    assert len(result) == 20  # noqa: PLR2004
    assert all(len(row) == 3 for row in result)  # noqa: PLR2004
    assert max(row[1] for row in result) > 0
    check.equal(result[0][0], pytest.approx(0.05, abs=1e-6))
    check.equal(result[0][1], pytest.approx(0.0, abs=1e-6))
    check.equal(result[14][1], pytest.approx(2.17714e-05, rel=1e-6))
    check.equal(result[14][2], pytest.approx(1.37695e-05, rel=1e-6))
    check.equal(result[-1][1], pytest.approx(0.0, abs=1e-6))


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


def test_traj_centering_for_maicos_mda(tip4p_trajectory, tmp_path):
    """Test that sequential MDA centering writes the expected trajectory."""
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.coordinates.XTC import XTCWriter

    from fairmd.lipids.analib.maicos import traj_centering_for_maicos_mda

    source, _ = tip4p_trajectory
    source.add_bonds([(0, 1), (0, 2), (0, 3)])
    topology = tmp_path / "topology.pdb"
    trajectory = tmp_path / "trajectory.xtc"
    source.atoms.write(topology)
    with XTCWriter(str(trajectory), n_atoms=source.atoms.n_atoms) as writer:
        for ts in source.trajectory:
            writer.write(source.atoms)

    universe = mda.Universe(topology, trajectory)
    output_dir = tmp_path / "centered"
    output_dir.mkdir()
    result = traj_centering_for_maicos_mda(
        universe,
        str(output_dir),
        "OW",
        eq_time=0,
        recompute=True,
    )

    assert os.path.isfile(result)
    assert result.endswith("whole.xtc")

    centered = mda.Universe(topology, result)
    assert centered.trajectory.n_frames == universe.trajectory.n_frames

    for ts in centered.trajectory:
        np.testing.assert_allclose(
            centered.select_atoms("name OW").center_of_mass()[2],
            ts.dimensions[2] / 2.0,
            atol=1e-4,
        )
