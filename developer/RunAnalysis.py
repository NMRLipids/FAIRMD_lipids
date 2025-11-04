"""
Execute pipeline of analysis methods that work globally on simulations.

.. note::
   This file is meant to be used by automated workflows.
"""

import logging
import subprocess
import sys

logging.basicConfig(
    format="%(asctime)s [%(levelname)s]: %(message)s",
    datefmt="%I:%M:%S %p",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)


def main() -> None:
    """Run analysis CLIs in sequence."""
    try:
        subprocess.run(["fmdl_match_experiments"], check=True)
        subprocess.run(["fmdl_evaluate_quality"], check=True)
        subprocess.run(["fmdl_make_ranking"], check=True)
    except subprocess.CalledProcessError:
        logger.exception("Run analysis failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
