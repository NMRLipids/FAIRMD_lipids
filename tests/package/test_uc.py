"""
Tests of `UniverseConstructor`, which builds MDAnalysis
universe and pre-download files.

Test folder: ToyData/Simulations.1

NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import copy
import os
import tempfile

import MDAnalysis as mda
import pytest
import pytest_check as check
# Global fixtures
# ----------------------------------------------------------------

pytestmark = [pytest.mark.sim1]


@pytest.fixture(scope="module")
def systems():
    from fairmd.lipids import FMDL_DATA_PATH
    from fairmd.lipids.core import initialize_databank

    if os.path.isfile(os.path.join(FMDL_DATA_PATH, ".notest")):
        pytest.exit("Test are corrupted. I see '.notest' file in the data folder.")
    s = initialize_databank()
    print(f"Loaded: {len(s)} systems")
    yield s


def test_uc_cleardata(systems, capsys):
    """Test if UC can clear predownloaded raw data"""
    from fairmd.lipids.api import UniverseConstructor

    s = systems.loc(566)
    uc = UniverseConstructor(s)

    # mocking that we downloaded struc-file
    tf = tempfile.NamedTemporaryFile(delete=False)
    tf.close()
    uc._paths["struc"] = tf.name

    uc.clear_mddata()
    check.is_in("struc-file", capsys.readouterr().out)
    check.is_none(uc.paths["struc"])
    check.is_false(os.path.isfile(tf.name))


def test_uc_build_wo_download(systems):
    """You must run `download` first even for predownloaded systems."""
    from fairmd.lipids.api import UniverseConstructor, UniverseConstructError

    for i in [243, 566]:
        s = systems.loc(i)
        uc = UniverseConstructor(s)
        with check.raises(UniverseConstructError):
            uc.build_universe()


def test_uc_skipdwnl_ifexists(systems):
    from fairmd.lipids import FMDL_SIMU_PATH
    from fairmd.lipids.api import UniverseConstructor

    s = systems.loc(566)
    uc = UniverseConstructor(s)
    for _sfx in ["tpr", "gro", "xtc"]:
        with open(os.path.join(FMDL_SIMU_PATH, s["path"], f"566.{_sfx}"), "w") as fd:
            fd.write(" ")
    uc.download_mddata()
    check.equal(os.stat(uc.paths["struc"]).st_size, 1)
    check.equal(os.stat(uc.paths["top"]).st_size, 1)
    check.equal(os.stat(uc.paths["traj"]).st_size, 1)
    # TEARDOWN
    uc.clear_mddata()


def test_download_and_clean_localhost(systems, capsys):
    """Check how is DOI=localhost system is processed."""
    from fairmd.lipids.api import UniverseConstructor

    s = systems.loc(243)
    uc = UniverseConstructor(s)
    # paths structure
    check.is_in("struc", uc.paths)
    check.is_in("top", uc.paths)
    check.is_in("traj", uc.paths)
    check.is_in("energy", uc.paths)

    uc.download_mddata()
    # md.gro is not downloaded
    check.is_true("md.gro" not in capsys.readouterr().err)
    # paths are filled in
    check.is_in(s["GRO"][0][0], uc.paths["struc"])
    check.is_true(os.path.isfile(uc.paths["struc"]))
    _store_fpath = uc.paths["struc"]
    uc.clear_mddata()
    check.is_none(uc.paths["struc"])
    check.is_true(os.path.isfile(_store_fpath))  # DOI-localhost systems are not really cleaned


@pytest.fixture(scope="function")
def fail_localhost_sys():
    from fairmd.lipids import FMDL_SIMU_PATH

    with tempfile.TemporaryDirectory(prefix=FMDL_SIMU_PATH + os.sep) as tmpd:
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
def test_fail_nonpredownld_localhost(fail_localhost_sys):
    from fairmd.lipids.api import UniverseConstructor

    uc = UniverseConstructor(fail_localhost_sys)
    uc.download_mddata()


def test_options_of_universe_construction(systems):
    from fairmd.lipids.api import UniverseConstructor
    from fairmd.lipids.core import System

    s0 = systems.loc(243)

    def get_nmols(s: System) -> int:
        _ans = 0
        for k, v in s["COMPOSITION"].items():
            if isinstance(v["COUNT"], list):
                _ans += sum(v["COUNT"])
            else:
                _ans += v["COUNT"]
        return _ans

    nmols = get_nmols(s0)

    # check default GMX downloader
    s = copy.deepcopy(s0)
    uc = UniverseConstructor(s)
    uc.download_mddata()  # GRO, XTC, TPR
    u = uc.build_universe()
    check.equal(u.atoms.molnums[-1] + 1, nmols)  # from-TPR universe has connectivity information
    check.greater(u.trajectory.n_frames, 1)
    uc.clear_mddata()

    # check what happends if we don't have TPR
    s = copy.deepcopy(s0)
    del s._store["TPR"]
    uc = UniverseConstructor(s)
    uc.download_mddata()  # now it will be "without" TPR
    u = uc.build_universe()
    with check.raises(mda.NoDataError):
        u.atoms.molnums  # Universe has no connectivity
    check.equal(u.atoms.n_atoms, s["NUMBER_OF_ATOMS"])  # but Universe is created
    check.greater(u.trajectory.n_frames, 1)
    uc.clear_mddata()


# 566 TPR GRO


@pytest.fixture(scope="function")
def uconstructor566(systems, capsys):
    from fairmd.lipids.api import UniverseConstructor

    s = systems.loc(566)
    uc = UniverseConstructor(s)
    check.equal(uc.system["ID"], 566)
    uc.download_mddata(skip_traj=True)
    assert "566.gro" in capsys.readouterr().out
    yield uc
    # TEARDOWN
    uc.clear_mddata()


def test_download_mddata_566_gro_tpr(uconstructor566):
    uc = uconstructor566
    s = uc.system

    # structure file is correct
    check.is_in(s["path"], uc.paths["struc"])
    check.is_in(s["GRO"][0][0], uc.paths["struc"])
    check.is_true(os.path.isfile(uc.paths["struc"]))

    # topology file is correct
    check.is_in(s["path"], uc.paths["top"])
    check.is_in(s["TPR"][0][0], uc.paths["top"])
    check.is_true(os.path.isfile(uc.paths["top"]))

    # traj is None for skip_traj=True
    check.is_none(uc.paths["traj"])


def test_universe_566_gro_tpr(uconstructor566):
    uc = uconstructor566
    s = uc.system
    u = uc.build_universe()
    check.is_true(isinstance(u, mda.Universe))
    check.equal(u.atoms.n_atoms, s["NUMBER_OF_ATOMS"])


def test_deal_with_corrupted_tpr(systems):
    from fairmd.lipids.api import UniverseConstructor
    from fairmd.lipids import FMDL_SIMU_PATH

    s = copy.deepcopy(systems.loc(86))
    with open(os.path.join(FMDL_SIMU_PATH, s["path"], f"86.tpr"), "w") as fd:
        fd.write("corrupted_data")

    # test without trajectory
    uc = UniverseConstructor(s)
    uc.download_mddata(skip_traj=True)
    u = uc.build_universe()  # build from GRO instead
    with check.raises(mda.NoDataError):
        u.atoms.molnums  # Universe has no connectivity
    check.equal(u.atoms.n_atoms, s["NUMBER_OF_ATOMS"])
    check.equal(u.trajectory.n_frames, 1)

    # test with the trajectory
    uc = UniverseConstructor(s)
    uc.download_mddata()
    u = uc.build_universe()
    with check.raises(mda.NoDataError):
        u.atoms.molnums  # Universe has no connectivity
    check.equal(u.atoms.n_atoms, s["NUMBER_OF_ATOMS"])
    check.greater(u.trajectory.n_frames, 1)
    uc.clear_mddata()


# 566 TPR XTC (no gro)


@pytest.fixture(scope="function")
def uconstructor566nogro(systems):
    from fairmd.lipids.api import UniverseConstructor

    s = copy.deepcopy(systems.loc(566))
    del s._store["GRO"]
    uc = UniverseConstructor(s)
    yield uc
    # TEARDOWN
    uc.clear_mddata()


def test_download_mddata_gromacs_tpr_xtc(uconstructor566nogro):
    uc = uconstructor566nogro
    s = uc.system
    uc.download_mddata()
    check.is_none(uc.paths["struc"])
