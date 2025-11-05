"""
Test module test_api.py

`test_api` contains tests of API those functions, which doesn't require making
MDAnalysis Universe and recomputing something. These functions just read README
files and precomputed JSON data.

Test data is stored in `./ToyData/Simulations.2`

-------------------------------------------------------------------------------
NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import os
import sys
import warnings
import pytest
import pytest_check as check

# run only on sim2 mocking data
pytestmark = [pytest.mark.sim2, pytest.mark.min]

# for vector comparisons with np.testing.assert_allclose
MAXRELERR_COMPARE_THRESHOLD = 1e-2
# testable constants
N_SYSTEMS_IN_TESTSET = 5


@pytest.fixture(scope="module")
def systems():
    """Fixture for loading the toy databank once per module."""
    from fairmd.lipids import FMDL_DATA_PATH
    from fairmd.lipids.core import initialize_databank

    if os.path.isfile(os.path.join(FMDL_DATA_PATH, ".notest")):
        pytest.exit("Test are corrupted. I see '.notest' file in the data folder.")
    s = initialize_databank()
    print(f"Loaded: {len(s)} systems")
    return s


def test_initialize_n(systems):
    """Check that we have 5 systems in the test set."""
    assert len(systems) == N_SYSTEMS_IN_TESTSET


def test_print_README(systems, capsys):
    from fairmd.lipids.core import print_README

    sys0 = systems[0]
    print_README(sys0)
    output: str = capsys.readouterr().out.rstrip()
    check.not_equal(output.find("DOI:"), -1)
    check.not_equal(output.find("TEMPERATURE:"), -1)
    print_README("example")
    output = capsys.readouterr().out.rstrip()
    check.is_in("Gromacs, Amber, NAMD", output)


@pytest.mark.parametrize(
    "systemid, result",
    [
        (281, 64.722),
        (566, 61.306),  # 1
        (787, 78.237),
        (243, 62.276),  # 2
        (86, 60.460),
    ],
)
def test_CalcAreaPerMolecule(systems, systemid, result):
    from fairmd.lipids.databankLibrary import CalcAreaPerMolecule

    sys0 = systems.loc(systemid)
    apm = CalcAreaPerMolecule(sys0)
    assert apm == pytest.approx(result, abs=6e-4)


@pytest.mark.parametrize(
    "systemid, result",
    [(281, 4142.234), (566, 3923.568), (787, 4694.191), (243, 2241.920), (86, 3869.417)],
)
def test_calcArea(systems, systemid, result):
    from fairmd.lipids.databankLibrary import calcArea

    sys0 = systems.loc(systemid)
    area = calcArea(sys0)
    assert area == pytest.approx(result, abs=6e-4)


@pytest.mark.parametrize("systemid, result", [(281, 128), (566, 128), (787, 120), (243, 72), (86, 128)])
def test_GetNLipids(systems, systemid, result):
    from fairmd.lipids.databankLibrary import GetNlipids

    sys0 = systems.loc(systemid)
    nlip = GetNlipids(sys0)
    assert nlip == result


@pytest.mark.parametrize(
    "systemid, result",
    [
        (281, [0.261, 0.415, 0.609]),
        (566, [0.260, 0.405, 0.597]),
        (243, [0.281, 0.423, 0.638]),
        (86, [0.264, 0.419, 0.623]),
    ],
)
def test_GetFormFactorMin(systems, systemid, result):
    import numpy as np

    from fairmd.lipids.databankLibrary import GetFormFactorMin

    sys0 = systems.loc(systemid)
    ffl = GetFormFactorMin(sys0)
    np.testing.assert_allclose(
        np.array(ffl[:3]),
        np.array(result),
        rtol=MAXRELERR_COMPARE_THRESHOLD,
        err_msg=(f"Problem in FFMIN comparison:\nComputed: {ffl[:3]!s} \nPre-computed: {result!s}"),
    )


@pytest.mark.parametrize("systemid, result", [(281, 31.5625), (566, 31.0), (787, 75.0), (243, 39.7778), (86, 27.75)])
def test_getHydrationLevel(systems, systemid, result):
    from fairmd.lipids.databankLibrary import getHydrationLevel

    sys0 = systems.loc(systemid)
    hl = getHydrationLevel(sys0)
    assert hl == pytest.approx(result, 1e-4)


@pytest.mark.parametrize(
    "systemid, lipid, result",
    [
        (281, ["POPC"], [1]),
        (566, ["POPC", "CHOL"], [0.9375, 0.0625]),
        (787, ["TOCL", "POPC", "POPE"], [0.25, 0.5, 0.25]),
        (243, ["DPPC"], [1]),
        (86, ["POPE"], [1]),
    ],
)
def test_calcLipidFraction(systems, systemid, lipid, result):
    from fairmd.lipids.databankLibrary import calcLipidFraction

    sys0 = systems.loc(systemid)
    assert calcLipidFraction(sys0, "SOPC") == 0  # absent lipid
    err = 0
    i = 0
    for i, lip in enumerate(lipid):
        err += (calcLipidFraction(sys0, lip) - result[i]) ** 2
    assert err == pytest.approx(0, rel=1e-4)


@pytest.mark.parametrize(
    "systemid, result",
    [
        (281, [-0.1610, -0.1217]),  # 1
        (566, [-0.1714, -0.1142]),  # 2
        (243, [-0.1764, -0.1784]),  # 3
        (86, [-0.1933, -0.1568]),
    ],
)
def test_averageOrderParameters(systems, systemid, result):
    from fairmd.lipids.databankLibrary import averageOrderParameters

    sys0 = systems.loc(systemid)
    sn1, sn2 = averageOrderParameters(sys0)
    assert (sn1 - result[0]) ** 2 + (sn2 - result[1]) ** 2 == pytest.approx(0, abs=1e-5)


# Tests behavior when averageOrderParameters cannot find calculated OP data


@pytest.mark.parametrize("systemid", [787])
def test_raises_averageOrderParameters(systems, systemid):
    from fairmd.lipids.databankLibrary import averageOrderParameters

    sys0 = systems.loc(systemid)
    with pytest.raises(FileNotFoundError) as exc_info:
        # sn1, sn2 = ..
        averageOrderParameters(sys0)
    assert "OrderParameters.json" in str(exc_info.value)


@pytest.mark.parametrize(
    "systemid, lipid, result",
    [
        (281, ["POPC/P"], ["M_G3P2_M"]),
        (566, ["POPC/P31", "CHOL/C1"], ["M_G3P2_M", "M_C1_M"]),
        (787, ["TOCL/P3", "POPC/P", "POPE/P"], ["M_G13P1_M", "M_G3P2_M", "M_G3P2_M"]),
        (243, ["DPPC/P8"], ["M_G3P2_M"]),
        (86, ["POPE/P8"], ["M_G3P2_M"]),
    ],
)
def test_getUniversalAtomName(systems, systemid, lipid, result):
    from fairmd.lipids.databankLibrary import getUniversalAtomName

    sys0 = systems.loc(systemid)
    i = 0
    for i, lipat in enumerate(lipid):
        lip, atom = tuple(lipat.split("/"))
        uname = getUniversalAtomName(sys0, atom, lip)
        check.equal(uname, result[i])


# Test fail-behavior of getUniversalAtomName


@pytest.mark.parametrize(
    "systemid, lipat, result",
    [(243, "DPPC/nonExisting", "Atom was not found"), (243, "nonExisting/P8", "was not found in the system")],
)
def test_bad_getUniversalAtomName(systems, systemid, lipat, result, capsys):
    from fairmd.lipids.databankLibrary import getUniversalAtomName

    sys0 = systems.loc(systemid)
    lip, atom = tuple(lipat.split("/"))
    uname = getUniversalAtomName(sys0, atom, lip)
    output = capsys.readouterr().err.rstrip()
    assert result in output
    assert uname is None


@pytest.mark.parametrize("systemid, lipid, result", [(243, "DPPC", "44ea5"), (787, "TOCL", "78629")])
def test_getAtoms(systems, systemid, lipid, result):
    from fairmd.lipids.databankLibrary import getAtoms

    sys0 = systems.loc(systemid)
    atoms = getAtoms(sys0, lipid).split()
    atoms = ",".join(sorted(atoms))
    import hashlib

    md5_hash = hashlib.md5()
    md5_hash.update(atoms.encode("ascii"))
    hx = md5_hash.hexdigest()[:5]
    assert hx == result


@pytest.mark.xfail(reason="Completely deprecated function", run=True, raises=NotImplementedError, strict=True)
@pytest.mark.parametrize("systemid, lipid, result", [(281, ["POPC"], [134])])
def test_raise_loadMappingFile(systems, systemid, lipid, result):
    from fairmd.lipids.databankLibrary import loadMappingFile

    sys0 = systems.loc(systemid)
    i = 0
    for i, lip in enumerate(lipid):
        mpf = loadMappingFile(sys0["COMPOSITION"][lip]["MAPPING"])
        check.equal(len(mpf), result[i])


@pytest.mark.parametrize(
    "systemid, lipid, result",
    [
        (281, ["POPC/P"], ["M_G3P2_M"]),
        (566, ["POPC/P31", "CHOL/C1"], ["M_G3P2_M", "M_C1_M"]),
        (787, ["TOCL/P3", "POPC/P", "POPE/P"], ["M_G13P1_M", "M_G3P2_M", "M_G3P2_M"]),
        (243, ["DPPC/P8"], ["M_G3P2_M"]),
        (86, ["POPE/P8"], ["M_G3P2_M"]),
    ],
)
def test_simulation2universal_atomnames(systems, systemid, lipid, result):
    from fairmd.lipids.databankLibrary import simulation2universal_atomnames

    sys0 = systems.loc(systemid)
    i = 0
    for i, lipat in enumerate(lipid):
        lip, atom = tuple(lipat.split("/"))
        sname = simulation2universal_atomnames(sys0, lip, result[i])
        check.equal(sname, atom)


@pytest.mark.parametrize(
    "systemid, lipat, result",
    [
        (243, "DPPC/nonExisting", "was not found from mappingDPPCberger.yaml"),
        (243, "nonExisting/M_G1_M", "was not found in the system!"),
    ],
)
def test_bad_simulation2universal_atomnames(systems, systemid, lipat, result, capsys):
    from fairmd.lipids.databankLibrary import simulation2universal_atomnames

    sys0 = systems.loc(systemid)
    lip, atom = tuple(lipat.split("/"))
    sname = simulation2universal_atomnames(sys0, lip, atom)
    output = capsys.readouterr().err.rstrip()
    assert result in output
    assert sname is None


@pytest.mark.parametrize(
    "systemid, result",
    [
        (281, "resname POPC"),
        (566, "resname CHL or resname OL or resname PA or resname PC"),
        (787, "resname POPC or resname POPE or resname TOCL2"),
        (243, "resname DPPC"),
        (86, "resname POPE"),
    ],
)
def test_getLipids(systems, systemid, result):
    from fairmd.lipids.databankLibrary import getLipids

    sys0 = systems.loc(systemid)
    gl = getLipids(sys0)
    assert gl == result


@pytest.fixture(scope="function")
def wipeth(systems):
    from fairmd.lipids import FMDL_SIMU_PATH

    # TD-FIXTURE FOR REMOVING THICKNESS AFTER TEST CALCULATIONS
    yield
    # TEARDOWN
    for sid in [243, 281]:
        sys0 = systems.loc(sid)
        fn = os.path.join(FMDL_SIMU_PATH, sys0["path"], "thickness.json")
        try:
            os.remove(fn)
        except FileNotFoundError:
            pass
        except Exception as e:
            print(f"An error occured during teardown: {e}", file=sys.stderr)
            raise


@pytest.mark.parametrize("systemid, result, thickres", [(281, 1, 4.19996), (243, 1, 4.25947)])
def test_analyze_th(systems, systemid, result, wipeth, thickres, logger):
    from fairmd.lipids import FMDL_SIMU_PATH
    from fairmd.lipids.analyze import computeTH

    sys0 = systems.loc(systemid)
    rc = computeTH(sys0, logger)
    assert rc == result
    fn = os.path.join(FMDL_SIMU_PATH, sys0["path"], "thickness.json")
    assert os.path.isfile(fn)
    with open(fn) as file:
        data = float(file.read().rstrip())
    assert data == pytest.approx(thickres, abs=1e-5)


@pytest.mark.parametrize("systemid, result", [(281, None), (243, None), (566, 4.2576), (787, None), (86, 4.1327)])
def test_GetThickness(systems, systemid, result):
    from fairmd.lipids.databankLibrary import GetThickness

    sys0 = systems.loc(systemid)
    th = GetThickness(sys0)
    if th is None and result is None:
        return
    assert float(th) == pytest.approx(float(result), abs=1e-4)


@pytest.mark.parametrize(
    "systemid, result",
    [
        (243, 0.7212884475213442),
        (86, 1.5018596337724872),
        (566, 1.1740608659926115),
    ],
)
def test_GetEquilibrationTimes(systems, systemid, result):
    from fairmd.lipids.databankLibrary import GetEquilibrationTimes
    from fairmd.lipids.settings.molecules import lipids_set

    sys0 = systems.loc(systemid)
    eq_times = GetEquilibrationTimes(sys0)
    print("\n========\n", eq_times, "\n=========\n")
    lips = list(set(sys0["COMPOSITION"].keys()).intersection(lipids_set.names))
    for lip in lips:
        if lip in ["CHOL", "DCHOL"]:
            continue  # for them eq_times are not computed
        check.is_in(lip, eq_times.keys())
        check.equal(result, eq_times[lip])


@pytest.mark.xfail(reason="EQtimes were not computed", run=True, raises=FileNotFoundError, strict=True)
def test_GetEquilibrationTimes_fail(systems):
    sys0 = systems.loc(787)
    from fairmd.lipids.databankLibrary import GetEquilibrationTimes

    GetEquilibrationTimes(sys0)


# Test that a valid JSON file is read correctly, left out full comparison of the OP-dictionary
@pytest.mark.parametrize(
    "systemid, lipid",
    [
        (281, "POPC"),
    ],
)
def test_GetOP_reads_valid_json(systems, systemid, lipid):
    from fairmd.lipids.databankLibrary import GetOP

    sys0 = systems.loc(systemid)
    resdic = GetOP(sys0)

    assert lipid in resdic


@pytest.mark.parametrize(
    "systemid, testmol ,result",
    [
        (787, "TOCL", None),
    ],
)
def test_GetOP_missing_file_warns(systems, systemid, testmol, result):
    from fairmd.lipids.databankLibrary import GetOP

    sys0 = systems.loc(systemid)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        resdic = GetOP(sys0)
        # Check the result is None for missing file
        assert testmol in resdic
        assert resdic[testmol] is result
        # Check a warning was raised containing the molecule name
        assert any(f"{testmol}OrderParameters.json not found" in str(wi.message) for wi in w)
