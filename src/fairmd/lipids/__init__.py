"""
FAIRMD Lipids project.

Here we define the main package variables and constants used throughout the FAIRMD
Lipids project. This includes *Package Information*, *Paths*, and *Return Codes*.
"""

import importlib.metadata
import importlib.util
import os
import re
import sys

from ._version import __version__

# Package Information
_mtd = importlib.metadata.metadata("fairmd-lipids")
__author_email__ = _mtd["Author-email"]
__author__ = re.sub(r" ?<.*>$", "", __author_email__)
__description__ = _mtd["Summary"]  # reads toml's `description` field
__license__ = _mtd["Licence-Expression"]
__url__ = _mtd.json["project_url"][2].split(", ")[1]  # url nr 3

# Global Paths


FMDL_DATA_PATH: str = os.environ.get(
    "FMDL_DATA_PATH",
    os.path.join(os.getcwd(), "BilayerData"),
)
""" Path to the Databank Data folder """

FMDL_SIMU_PATH: str = os.environ.get(
    "FMDL_SIMU_PATH",
    os.path.join(FMDL_DATA_PATH, "Simulations"),
)
""" Path to the project simulations folder"""

FMDL_MOL_PATH: str = os.path.join(FMDL_DATA_PATH, "Molecules")
""" Path to the project molecules folder """

FMDL_EXP_PATH: str = os.path.join(FMDL_DATA_PATH, "experiments")
""" Path to the project experiments folder """

# Universal return codes

RCODE_SKIPPED: int = 0
""" Success code 0: calculation skipped """

RCODE_COMPUTED: int = 1
""" Success code 1: calculation successful """

RCODE_ERROR: int = 2
""" Success code 2: calculation failed """


print(
    f"fairmd-lipids {__version__} by {__author__} - {__license__}",
)

if os.path.isdir(FMDL_DATA_PATH):

    def raise_if_subpath_of_dblspec(p: str) -> None:
        """Raise if p is a subpath of the package path."""
        pkgpth = importlib.util.find_spec("fairmd.lipids").submodule_search_locations[0]
        ndp = os.path.abspath(p)
        if os.path.commonpath([pkgpth, ndp]) == pkgpth:
            msg = (
                "We forbid using ToyData in the distribution package directly."
                "Please run `fmdl_initialize_data toy` instead"
            )
            raise RuntimeError(msg)

    for p in [FMDL_DATA_PATH, FMDL_EXP_PATH, FMDL_MOL_PATH, FMDL_SIMU_PATH]:
        raise_if_subpath_of_dblspec(p)

    from fairmd.lipids.settings import molecules

    _ = len(molecules.lipids_set)
    print(
        f"FAIRMD Lipids is initialized from the folder: {FMDL_DATA_PATH}\n"
        "---------------------------------------------------------------",
    )
elif "fmdl_initialize_data" in sys.argv[0]:
    # fmdl_initialize_data is used to create databank directories
    # so we should not complain that directories don't exist
    pass
else:
    msg = f"""
Error: no data folder {FMDL_DATA_PATH}.
If Data folder was not created, please create it by using
 $ fmdl_initialize_data.py toy
          OR
 $ fmdl_initialize_data.py stable
and then specify by FMDL_DATA_PATH environment variable."""
    raise RuntimeError(msg)


__all__ = [
    "FMDL_DATA_PATH",
    "FMDL_EXP_PATH",
    "FMDL_MOL_PATH",
    "FMDL_SIMU_PATH",
    "RCODE_COMPUTED",
    "RCODE_ERROR",
    "RCODE_SKIPPED",
    "__author__",
    "__author_email__",
    "__description__",
    "__license__",
    "__url__",
    "__version__",
]
