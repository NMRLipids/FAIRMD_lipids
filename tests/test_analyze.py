"""
`test_analyze` tests functions performing databank recomputing. It starts mostly from
downloading everything from testing repository.

Test data is stored in `./ToyData/Simulations.1`

NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import glob
import logging
import os

import pytest
import pytest_check as check

# Set up logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# Create formatter and add it to the handlers
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)

# Add the handler to the logger
logger.addHandler(ch)

# Global loaders
# ----------------------------------------------------------------
# we substitute simulation path to test databank analysis scripts
# we load all trajectories from specific zenodo repository
# https://doi.org/10.5281/zenodo.11614468
# On teardown stage we clear the folders.

# run only on sim2 mocking data
pytestmark = pytest.mark.sim1


@pytest.fixture(scope="module")
def systems():
    import fairmd.lipids
    from fairmd.lipids.core import initialize_databank

    logger.info("Starting systems fixture initialization")
    if os.path.isfile(os.path.join(lipids.NMLDB_DATA_PATH, ".notest")):
        pytest.exit("Test are corrupted. I see '.notest' file in the data folder.")
    s = initialize_databank()
    logger.info(f"Loaded: {len(s)} systems")
    yield s
    # TEARDOWN SYSTEMS
    logger.info("Wiping temporary calculation data.")
    for _s in s:

        def gbGen(x):
            return glob.glob(os.path.join(lipids.NMLDB_SIMU_PATH, _s["path"], x))

        clear_list = []
        if "NMLDB_TEST_NOWIPE" in os.environ:
            # leave JSONs
            pass
        else:
            clear_list = ["*.json"]
        clear_list += [
            "*.dat",
            "conf.gro",
            "frame0.gro",
            "*.dat",
            "*.buildH",
            ".*",
            "#*",
            "*.def",
            "whole.xtc",
            "centered.xtc",
        ]
        for pat in clear_list:
            for f in gbGen(pat):
                logger.debug(f"Removing file: {f}")
                os.remove(f)


@pytest.fixture(scope="module")
def systemLoadTraj(systems):
    import fairmd.lipids
    from fairmd.lipids.databankLibrary import system2MDanalysisUniverse

    logger.info("Downloading trajectory data.")
    for s in systems:
        logger.debug(f"Processing system: {s.get('id', 'unknown')}")
        u = system2MDanalysisUniverse(s)
        del u
    yield
    # TEARDOWN SYSTEM-LOADING
    if "NMLDB_TEST_NOWIPE" in os.environ:
        logger.info("DBG: Skipping trajectory data wipe.")
    else:
        print("DBG: Wiping trajectory data.")
        for s in systems:
            logger.debug("Skipping localhost system")
            if s.get("DOI") == "localhost":
                continue  # don't remove in this case!
            for wp in ["GRO", "TPR", "TRJ"]:
                try:
                    file_ = s[wp][0][0]
                except (KeyError, TypeError):
                    logger.warning(f"No {wp} file found for system")
                    continue
                file_ = os.path.join(lipids.NMLDB_SIMU_PATH, s["path"], file_)
                logger.debug(f"Checking file: {file_}")
                if os.path.exists(file_):
                    os.remove(file_)


# Test functions block.
# ----------------------------------------------------------------
# Every test function is parametrized with system ID to make clear reporting
# about which system actually fails on a test function.
MAXRELERR_COMPARE_THRESHOLD = 1e-2


def compareJSONsBtwSD(jsfn: str):
    """Function to call from the test to compare JSON to the precomputed one.
    Raises standard assertion to fail the test.

    Args:
        jsfn (_type_): relative path from simulation folder
    """
    import json
    import numpy as np
    from fairmd.lipids import NMLDB_DATA_PATH

    logger.info(f"Comparing JSON files for: {jsfn}")

    _p1 = os.path.join(NMLDB_DATA_PATH, "Simulations.1")
    _p2 = os.path.join(NMLDB_DATA_PATH, "Simulations.2")
    jsf1 = os.path.join(_p1, jsfn)
    jsf2 = os.path.join(_p2, jsfn)

    with open(jsf1) as f:
        j1 = json.load(f)
    with open(jsf2) as f:
        j2 = json.load(f)

    assert type(j1) is type(j2)

    if type(j1) is list:
        assert len(j1) == len(j2), f"Problem in {jsfn} comparison: lists has different lengths!"
        # to process further list as dict
        j1 = dict(enumerate(j1))
        j2 = dict(enumerate(j2))

    if type(j1) is dict:
        assert set(j1.keys()) == set(j2.keys())
        for k1 in j1:
            np.testing.assert_allclose(
                np.array(j1[k1]),
                np.array(j2[k1]),
                rtol=MAXRELERR_COMPARE_THRESHOLD,
                err_msg=(f"Problem in {jsfn} comparison: {k1} line.\nComputed: {j1[k1]!s} \nPre-computed: {j2[k1]!s}"),
            )

    logger.info(f"Data {jsfn} was compared against precomputed!")


@pytest.mark.parametrize("systemid", [86, 243, 281, 566, 787])
def test_analyze_apl(systems, systemLoadTraj, systemid, logger):
    import fairmd.lipids
    from fairmd.lipids.analyze import computeAPL

    s = systems.loc(systemid)

    rCode = computeAPL(s, logger)
    assert rCode == lipids.analyze.RCODE_COMPUTED

    cFile = os.path.join(lipids.NMLDB_SIMU_PATH, s["path"], "apl.json")
    assert os.path.isfile(cFile)
    assert os.path.getsize(cFile) > 1e3
    compareJSONsBtwSD(os.path.relpath(cFile, lipids.NMLDB_SIMU_PATH))


@pytest.mark.parametrize(
    "systemid, rcodex",
    [(281, 1), (566, 1), (787, 2), (243, 1), (86, 1)],
)
def test_analyze_op(systems, systemLoadTraj, systemid, rcodex, logger):
    import fairmd.lipids
    from fairmd.lipids.analyze import computeOP
    from fairmd.lipids.settings.molecules import lipids_set

    s = systems.loc(systemid)
    rCode = computeOP(s, logger)
    assert rCode == rcodex, f"rCode = {rCode} is not the same as rcodex = {rcodex}"
    if rcodex == lipids.RCODE_ERROR:
        return

    for lip in set(s["COMPOSITION"].keys()).intersection(set(lipids_set.names)):
        cFile = os.path.join(
            lipids.NMLDB_SIMU_PATH,
            s["path"],
            lip + "OrderParameters.json",
        )
        assert os.path.isfile(cFile), f"File {cFile} wasn't created for {lip}!"
        assert os.path.getsize(cFile) > 1e3, f"File {cFile} for {lip} is less than 1K!"
        compareJSONsBtwSD(os.path.relpath(cFile, lipids.NMLDB_SIMU_PATH))


@pytest.mark.parametrize("systemid, rcodex", [(281, 1), (566, 1), (787, 2), (243, 1), (86, 1)])
def test_analyze_maicos(systems, systemLoadTraj, systemid, rcodex, logger):
    import fairmd.lipids
    from fairmd.lipids.analyze import computeMAICOS

    s = systems.loc(systemid)
    rCode = computeMAICOS(s, logger)
    assert rCode == rcodex, f"rCode = {rCode} is not the same as rcodex = {rcodex}"
    if rcodex == lipids.RCODE_ERROR:
        return

    for fn in ["WaterDensity.json", "LipidDensity.json", "TotalDensity.json", "FormFactor.json"]:
        cFile = os.path.join(lipids.NMLDB_SIMU_PATH, s["path"], fn)
        check.is_true(os.path.isfile(cFile))
        check.is_true(os.path.getsize(cFile) > 1e3)
        try:
            compareJSONsBtwSD(
                os.path.relpath(cFile, lipids.NMLDB_SIMU_PATH),
            )
        except AssertionError as e:
            check.is_true(False, msg=str(e))
    # TODO: check other MAICoS outputs for at least existence


@pytest.mark.parametrize(
    "systemid, rcodex",
    [(281, 1), (566, 1), (787, 2), (243, 1), (86, 1)],
)
def test_analyze_nmrpca(systems, systemLoadTraj, systemid, rcodex, logger):
    import fairmd.lipids
    from fairmd.lipids.analyze import computeNMRPCA

    s = systems.loc(systemid)
    rCode = computeNMRPCA(s, logger)
    assert rCode == rcodex

    if rCode == lipids.RCODE_ERROR:
        return

    cFile = os.path.join(lipids.NMLDB_SIMU_PATH, s["path"], "eq_times.json")
    assert os.path.isfile(cFile)
    assert os.path.getsize(cFile) > 10
    compareJSONsBtwSD(os.path.relpath(cFile, lipids.NMLDB_SIMU_PATH))
