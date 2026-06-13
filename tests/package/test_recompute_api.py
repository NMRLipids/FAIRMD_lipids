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
    from fairmd.lipids.api import mda_read_trj_tilt_angles, UniverseConstructor

    s = systems.loc(systemid)
    uc = UniverseConstructor(s)
    uc.download_mddata()
    u = uc.build_universe()
    pats = u.select_atoms(s.content[lipid].uan2selection("M_G3P2_M", lipid))
    nats = u.select_atoms(s.content[lipid].uan2selection("M_G3N6_M", lipid))
    a1 = pats.atoms.names[0]
    a2 = nats.atoms.names[0]

    a, b, c, d = mda_read_trj_tilt_angles(lipid, a1, a2, u)
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
