"""
`test_recompute_api.py` tests those parts of API, which require building MDAnalysis
universe and from-trajectories computations.

Test folder: ToyData/Simulations.1

NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import glob
import os
from contextlib import nullcontext as does_not_raise
from tempfile import TemporaryDirectory

import pytest

# Global fixtures
# ----------------------------------------------------------------

pytestmark = [pytest.mark.sim1]


@pytest.fixture(scope="module")
def systems():
    from fairmd.lipids import FMDL_DATA_PATH, FMDL_SIMU_PATH
    from fairmd.lipids.core import initialize_databank

    if os.path.isfile(os.path.join(FMDL_DATA_PATH, ".notest")):
        pytest.exit("Test are corrupted. I see '.notest' file in the data folder.")
    s = initialize_databank()
    print(f"Loaded: {len(s)} systems")
    yield s
    # TEARDOWN SYSTEMS
    print("DBG: Wiping temporary calculation data.")
    for _sid in [787]:
        _s = s.loc(_sid)

        def gbGen(x):
            return glob.glob(os.path.join(FMDL_SIMU_PATH, _s["path"], x))

        clearList = ["*.xtc", "*.gro"]
        for pat in clearList:
            for f in gbGen(pat):
                os.remove(f)


# Test functions block.
# ----------------------------------------------------------------
# Every test function is parametrized with system ID to make clear reporting
# about which system actually fails in a test function.


@pytest.mark.parametrize(
    "systemid, natoms, nframes",
    [
        (243, 12208, 121),  # with TPR, united-atom, local
        (787, 46290, 1141),  # with GRO only, aa, network
    ],
)
def test_system2MDAnalysisUniverse(systems, systemid, natoms, nframes):
    from fairmd.lipids.databankLibrary import system2MDanalysisUniverse

    s = systems.loc(systemid)
    u = system2MDanalysisUniverse(s)
    assert len(u.atoms) == natoms
    assert u.trajectory.n_frames == nframes
    with does_not_raise() as _:
        # check that it doesn't fail iterating over the trajectory
        for _ in u.trajectory:
            pass


@pytest.fixture(scope="function")
def failSys():
    from fairmd.lipids import FMDL_SIMU_PATH

    with TemporaryDirectory(prefix=FMDL_SIMU_PATH + os.sep) as tmpd:
        s = {
            "DOI": "localhost",
            "GRO": [["md.gro"]],
            "TRJ": [["md.trr"]],
            "path": os.path.relpath(FMDL_SIMU_PATH, tmpd),
            "SOFTWARE": "gromacs",
        }
        yield s
    # TEARDOWN


@pytest.mark.xfail(reason="Localhost with non-downloaded files", raises=FileNotFoundError)
def test_fail1_system2MDAnalysisUniverse(failSys):
    from fairmd.lipids.databankLibrary import system2MDanalysisUniverse

    _ = system2MDanalysisUniverse(failSys)


def hashFV(x):
    import numpy as np

    a = np.array(x)
    a = np.around(a, 4)
    a *= 1e4
    a = np.array(a, dtype="int32")
    a = tuple(a.tolist())
    return hash(a)


@pytest.mark.parametrize(
    "systemid, lipid, fvhash",
    [
        (243, "DPPC", -5227956720741036084),  # with TPR, united-atom, local
        (787, "POPC", 4799549858726566566),  # with GRO only, aa, network
    ],
)
def test_PJangle(systems, systemid, lipid, fvhash):
    from fairmd.lipids.databankLibrary import (
        read_trj_PN_angles,
        simulation2universal_atomnames,
        system2MDanalysisUniverse,
    )

    s = systems.loc(systemid)
    u = system2MDanalysisUniverse(s)
    a1 = simulation2universal_atomnames(s, lipid, "M_G3P2_M")
    a2 = simulation2universal_atomnames(s, lipid, "M_G3N6_M")

    a, b, c, d = read_trj_PN_angles(lipid, a1, a2, u)
    # time-molecule arrays
    assert len(a) == sum(s["COMPOSITION"][lipid]["COUNT"])
    # time-averaged list
    assert len(b) == sum(s["COMPOSITION"][lipid]["COUNT"])
    # comparing per-molecule hash. Suppose that's enough
    assert hashFV(b) == fvhash
    # overall mean
    assert c > 0
    # overall std
    assert d > 0
