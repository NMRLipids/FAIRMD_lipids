"""
Validate integrity of NMRlipids Data.

This script verifies that the databank can be initialized and that lipid
and molecule sets can be constructed from the data files.
"""

import logging
import sys

logging.basicConfig(
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%I:%M:%S %p",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)


def check_integrity() -> None:
    """Run imports of lipid_set and molecule_set"""
    from fairmd.lipids.core import initialize_databank
    from fairmd.lipids.molecules import lipids_set, molecules_set  # noqa: F401

    initialize_databank()


if __name__ == "__main__":
    try:
        check_integrity()
        logger.info("fairmd-lipids integrity check passed")
    except Exception:
        logger.exception("fairmd-lipids integrity check failed")
        sys.exit(1)
