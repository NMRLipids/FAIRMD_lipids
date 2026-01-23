import os
import pytest
import logging
import sys

# Pytest HOOKS
# -------------------------------------------------------------------


def pytest_addoption(parser):
    parser.addoption(
        "--cmdopt",
        action="store",
        default="sim1",
        help="Two test groups: sim1|sim2|nodata",
    )


# Pytest GLOBAL FIXTURES
# -------------------------------------------------------------------

# Dictionary for pytest marker to simulation folder
SIM_MAP = {
    "sim1": "Simulations.1",
    "sim2": "Simulations.2",
    "adddata": "Simulations.AddData",
    "rdkit": "Simulations.2",
    "nodata": None,
}


@pytest.fixture(autouse=True, scope="module")
def header_module_scope(request):
    """Set env vars depending on data required, remove fairmd-lipids on teardown."""
    # Find pytest marker:
    sim_key = None
    for key in ("sim1", "sim2", "adddata", "nodata"):
        if request.node.get_closest_marker(key):
            sim_key = key
            break

    # Default to nodata:
    if sim_key is None:
        cmdopt = request.config.getoption("--cmdopt")
        sim_key = cmdopt if cmdopt in SIM_MAP else "nodata"

    data_root = os.path.join(os.path.dirname(__file__), "ToyData")
    os.environ["FMDL_DATA_PATH"] = data_root
    sim_dir = SIM_MAP[sim_key]
    if sim_dir:
        os.environ["FMDL_SIMU_PATH"] = os.path.join(data_root, sim_dir)
    else:
        os.environ.pop("FMDL_SIMU_PATH", None)

    print("DBG env -> FMDL_DATA_PATH:", os.getenv("FMDL_DATA_PATH"))
    print("DBG env -> FMDL_SIMU_PATH:", os.getenv("FMDL_SIMU_PATH"))

    yield
    # Teardown:
    remove_databank_import()
    print("DBG: Mocking completed")


@pytest.fixture(scope="module")
def logger():
    logger = logging.getLogger("test_logger")
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)

    if not logger.handlers:  # Avoid adding multiple handlers during pytest runs
        logger.addHandler(handler)

    yield logger

    # TEARDOWN: clean up handlers after use
    for handler in logger.handlers:
        logger.removeHandler(handler)


def remove_databank_import() -> None:
    """Delete fairmd.lipids module from sys.modules, resets future imports"""
    for name in list(sys.modules):
        if name.startswith("fairmd.lipids"):
            del sys.modules[name]
