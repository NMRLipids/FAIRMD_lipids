"""
Test module for ranking table creation.

Test data is stored in `./ToyData/Simulations.2`

-------------------------------------------------------------------------------
NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import os
import shutil
import pytest
import pytest_check as check

# run only on sim2 mocking data
pytestmark = [pytest.mark.sim2, pytest.mark.min]


def test_make_ranking():
    """Test the make_ranking function."""
    from fairmd.lipids.bin.make_ranking import make_ranking
    from fairmd.lipids import FMDL_DATA_PATH

    rnkpath = os.path.join(FMDL_DATA_PATH, "Ranking")
    try:
        os.makedirs(rnkpath)

        make_ranking()

        check.is_true(os.path.isfile(os.path.join(rnkpath, "FF_ranking.csv")), "FF_ranking.csv not created")
        check.is_true(os.path.isfile(os.path.join(rnkpath, "OP_ranking.csv")), "OP_ranking.csv not created")
        check.is_true(os.path.isfile(os.path.join(rnkpath, "POPC_ranking.csv")), "OP_ranking.csv not created")
    finally:
        shutil.rmtree(rnkpath)
