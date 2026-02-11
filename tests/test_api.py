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

import copy
import os
import sys
import warnings
import numpy as np
import pytest
import pytest_check as check

# run only on sim2 mocking data
pytestmark = [pytest.mark.sim2, pytest.mark.min]

# for vector comparisons with np.testing.assert_allclose
MAXRELERR_COMPARE_THRESHOLD = 1e-2
# testable constants
N_SYSTEMS_IN_TESTSET = 6


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


def test_content(systems):
    from fairmd.lipids.molecules import Lipid, Molecule

    for s in systems:
        check.is_instance(s.content, dict)
        check.is_instance(s.lipids, dict)
        check.is_instance(s.n_lipids, int)
        check.equal(len(s.content), len(s["COMPOSITION"]))  # incl. water
        for k, v in s.content.items():
            check.is_in(k, s["COMPOSITION"])
            check.is_instance(v, Molecule)
        for k, v in s.lipids.items():
            check.is_in(k, s["COMPOSITION"])
            check.is_instance(v, Lipid)


def test_copy_system(systems):
    s1 = systems.loc(566)
    s2 = copy.deepcopy(s1)
    s2._store["GRO"] = None
    assert s2["GRO"] is None
    assert s1["GRO"] is not None


def test_mda_gen_selection_mols(systems):
    from fairmd.lipids.api import mda_gen_selection_mols
    from fairmd.lipids.molecules import lipids_set

    sys566 = systems.loc(566)
    alllips = mda_gen_selection_mols(sys566)
    check.equal(alllips, "resname CHL or resname OL or resname PA or resname PC")
    pconly = mda_gen_selection_mols(sys566, molecules=[lipids_set.get("POPC")])
    check.equal(pconly, "resname OL or resname PA or resname PC")

    sys787 = systems.loc(787)
    alllips = mda_gen_selection_mols(sys787)
    check.equal(alllips, "resname POPC or resname POPE or resname TOCL2")
    pconly = mda_gen_selection_mols(sys787, molecules=[lipids_set.get("POPC")])
    check.equal(pconly, "resname POPC")


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
def test_get_mean_apl(systems, systemid, result):
    from fairmd.lipids.api import get_mean_ApL

    sys0 = systems.loc(systemid)
    apm = get_mean_ApL(sys0)
    assert apm == pytest.approx(result, abs=6e-4)


@pytest.mark.parametrize(
    "systemid, result",
    [(281, 4142.234), (566, 3923.568), (787, 4694.191), (243, 2241.920), (86, 3869.417)],
)
def test_get_total_area(systems, systemid, result):
    from fairmd.lipids.api import get_total_area

    sys0 = systems.loc(systemid)
    area = get_total_area(sys0)
    assert area == pytest.approx(result, abs=6e-4)


@pytest.mark.parametrize("systemid, result", [(281, 128), (566, 128), (787, 120), (243, 72), (86, 128)])
def test_n_lipids(systems, systemid, result):
    sys0 = systems.loc(systemid)
    nlip = sys0.n_lipids
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

    from fairmd.lipids.api import get_FF
    from fairmd.lipids.analib.formfactor import get_mins_from_ffdata

    sys0 = systems.loc(systemid)
    ff_data = get_FF(sys0)
    ffl = get_mins_from_ffdata(ff_data)
    np.testing.assert_allclose(
        np.array(ffl[:3]),
        np.array(result),
        rtol=MAXRELERR_COMPARE_THRESHOLD,
        err_msg=(f"Problem in FFMIN comparison:\nComputed: {ffl[:3]!s} \nPre-computed: {result!s}"),
    )


@pytest.mark.parametrize("systemid, result", [(281, 31.5625), (566, 31.0), (787, 75.0), (243, 39.7778), (86, 27.75)])
def test_get_hydration(systems, systemid, result):
    sys0 = systems.loc(systemid)
    hl = sys0.get_hydration()
    check.almost_equal(hl, sys0.get_hydration(basis="number"))  # number is default
    check.almost_equal(hl, result, abs=1e-4)


@pytest.mark.parametrize(
    "systemid, lipid, result_molar, result_mass",
    [
        (281, ["POPC"], [1], [1]),
        (566, ["POPC", "CHOL"], [0.9375, 0.0625], [0.9672, 0.0328]),
        (787, ["TOCL", "POPC", "POPE"], [0.25, 0.5, 0.25], [0.3945, 0.4113, 0.1941]),
        (243, ["DPPC"], [1], [1]),
        (86, ["POPE"], [1], [1]),
    ],
)
def test_membrane_composition(systems, systemid, lipid, result_molar, result_mass):
    sys0 = systems.loc(systemid)
    molar_fractions = sys0.membrane_composition(basis="molar")
    mass_fractions = sys0.membrane_composition(basis="mass")
    with check.raises(ValueError):
        _ = sys0.membrane_composition(basis="invalid_option")
    with check.raises(KeyError):
        _ = molar_fractions["SOPC"]
    err = 0
    for i, lip in enumerate(lipid):
        err += (molar_fractions[lip] - result_molar[i]) ** 2
    check.almost_equal(err, 0, abs=1e-5)
    err = 0
    for i, lip in enumerate(lipid):
        err += (mass_fractions[lip] - result_mass[i]) ** 2
    check.almost_equal(err, 0, abs=1e-5)


@pytest.mark.parametrize(
    "systemid, result_molar",
    [
        (243, {"SOD": 0.155, "CLA": 0.155}),
        (787, {"SOD": 0.37}),
        (281, {}),
    ],
)
def test_solution_composition(systems, systemid, result_molar):
    sys0 = systems.loc(systemid)
    molar_fractions = sys0.solution_composition(basis="molar")
    with check.raises(ValueError):
        _ = sys0.solution_composition(basis="invalid_option")
    check.equal(molar_fractions.keys(), result_molar.keys())
    for k, v in result_molar.items():
        check.almost_equal(molar_fractions[k], v, abs=1e-3)


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
    from fairmd.lipids.api import mda_gen_selection_mols

    sys0 = systems.loc(systemid)
    gl = mda_gen_selection_mols(sys0)
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


@pytest.mark.parametrize(
    "systemid, result",
    [
        (281, None),
        (243, None),
        (787, None),  # exceptions are expected
        (566, 4.2576),
        (86, 4.1327),
    ],
)
def test_get_thickness(systems, systemid, result):
    from fairmd.lipids.api import get_thickness
    import numbers

    sys0 = systems.loc(systemid)
    if result is None:
        with pytest.raises(FileNotFoundError):
            get_thickness(sys0)
        return
    th = get_thickness(sys0)
    check.is_true(isinstance(th, numbers.Real))
    check.almost_equal(th, result, abs=1e-4)


@pytest.mark.parametrize(
    "systemid, result",
    [
        (243, 0.7212884475213442),
        (86, 1.5018596337724872),
        (566, 1.1740608659926115),
        (787, None),
    ],
)
def test_get_eqtimes(systems, systemid, result):
    from fairmd.lipids.api import get_eqtimes
    from fairmd.lipids.molecules import lipids_set

    sys0 = systems.loc(systemid)

    if result is None:
        with pytest.raises(FileNotFoundError):
            get_eqtimes(sys0)
        return

    eq_times = get_eqtimes(sys0)
    print("\n========\n", eq_times, "\n=========\n")
    lips = list(set(sys0["COMPOSITION"].keys()).intersection(lipids_set.names))
    for lip in lips:
        if lip in ["CHOL", "DCHOL"]:
            continue  # for them eq_times are not computed
        check.is_in(lip, eq_times.keys())
        check.equal(result, eq_times[lip])


# Test that a valid JSON file is read correctly, left out full comparison of the OP-dictionary
@pytest.mark.parametrize(
    "systemid, lipid",
    [
        (281, "POPC"),
    ],
)
def test_get_OP_reads_valid_json(systems, systemid, lipid):
    from fairmd.lipids.api import get_OP

    sys0 = systems.loc(systemid)
    resdic = get_OP(sys0)

    assert lipid in resdic
    np.testing.assert_allclose(resdic[lipid]["M_G1_M M_G1H1_M"], np.array([-0.169826, 0.0268957, 0.00238661]))


@pytest.mark.parametrize(
    "systemid",
    [
        281,
        243,
    ],
)
def test_get_FF_valid(systems, systemid):
    from fairmd.lipids.api import get_FF

    sys0 = systems.loc(systemid)
    ff_data: np.ndarray = get_FF(sys0)

    assert isinstance(ff_data, np.ndarray)
    check.equal(ff_data.shape[1], 3)
    check.greater(ff_data.shape[0], 10)
    check.is_true((ff_data[:, 0] < 1.01).all(), "First column is qq")
    check.greater(np.mean(ff_data[:, 1] / ff_data[:, 2]), 1, "Error should be smaller than values")


@pytest.mark.parametrize(
    "systemid",
    [
        787,
    ],
)
def test_get_FF_missing_file(systems, systemid):
    from fairmd.lipids.api import get_FF

    sys0 = systems.loc(systemid)

    with check.raises(FileNotFoundError) as e:
        ff_data = get_FF(sys0)
    check.is_in("form-factor data is missing", str(e.value))


@pytest.mark.parametrize(
    "systemid, testmol ,result",
    [
        (787, "TOCL", None),
    ],
)
def test_GetOP_missing_file_warns(systems, systemid, testmol, result):
    from fairmd.lipids.api import get_OP

    sys0 = systems.loc(systemid)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        resdic = get_OP(sys0)
        # Check the result is None for missing file
        assert testmol in resdic
        assert resdic[testmol] is result
        # Check a warning was raised containing the molecule name
        assert any(f"{testmol}OrderParameters.json not found" in str(wi.message) for wi in w)


def test_run_analysis_interface():
    from fairmd.lipids import RCODE_COMPUTED, RCODE_ERROR, RCODE_SKIPPED
    from fairmd.lipids.utils import run_analysis

    # create logger for testing
    import logging
    from io import StringIO

    log_stream = StringIO()
    _logger = logging.getLogger("test_run_analysis")
    _logger.setLevel(logging.INFO)
    _logger.addHandler(logging.StreamHandler(log_stream))

    def dummy_method(system, logger):
        logger.info(f"Dummy method called for system ID {system['ID']}")
        return RCODE_COMPUTED

    run_analysis(
        method=dummy_method,
        logger=_logger,
        id_range=(None, None),
    )

    check.is_in(f"COMPUTED: {N_SYSTEMS_IN_TESTSET}", log_stream.getvalue())
    check.is_in("SKIPPED: 0", log_stream.getvalue())

    run_analysis(
        method=dummy_method,
        logger=_logger,
        id_list=[86, 281, 243],
    )

    check.is_in("COMPUTED: 3", log_stream.getvalue())


def test_get_quality(systems):
    from fairmd.lipids.api import get_quality

    sys0 = systems.loc(281)
    # default quality is TOTAL|BOTH
    # "boths" experiments work only for "total"
    q1: float = 0
    q2: float = 0
    with check.raises(NotImplementedError):
        q1 = get_quality(sys0)
    with check.raises(NotImplementedError):
        q2 = get_quality(sys0, part="total", experiment="both")
    # UNCOMMENT WHEN IMPLEMENTED
    check.equal(q1, q2, "Default run should be TOTAL|BOTH")
    # for parts "both" is not allowed
    with check.raises(ValueError):
        get_quality(sys0, part="headgroup", experiment="both")
    with check.raises(ValueError):
        get_quality(sys0, lipid="POPC", experiment="both")

    # total quality can work for each experiment
    q1 = get_quality(sys0, part="total", experiment="FF")
    q2 = get_quality(sys0, experiment="FF")  # part=total is default
    check.equal(q1, q2, "Default `part` must be TOTAL")
    check.equal(q1, 0.8, "FormFactor quality value is improper!")
    q = get_quality(sys0, part="total", experiment="OP")
    check.equal(q, 0.41)

    # FF quality for parts is not allowed
    with check.raises(ValueError):
        get_quality(sys0, part="headgroup", experiment="FF")
    with check.raises(ValueError):
        get_quality(sys0, lipid="POPC", experiment="FF")
    # NMR quality can work for parts
    q = get_quality(sys0, part="headgroup", experiment="OP")
    check.equal(q, 0.5)
    q = get_quality(sys0, lipid="POPC", experiment="OP")
    check.equal(q, 0.4)
    q = get_quality(sys0, part="tails", lipid="POPC", experiment="OP")
    check.equal(q, 0.25)
    # invalid part
    with check.raises(ValueError):
        get_quality(sys0, part="invalidpart", experiment="OP")
    # invalid lipid
    with check.raises(ValueError):
        get_quality(sys0, lipid="INVALIDLIPID", experiment="OP")
    sys0 = systems.loc(566)
    q = get_quality(sys0, part="tails", lipid="POPC", experiment="OP")
    check.is_nan(q, "Absent quality must mean nan")
