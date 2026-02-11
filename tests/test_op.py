"""
Test Order Parameter handling functionality.

Currently focused on testing opconvertor realted code because it is important for
OP data visualization.

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


class TestBuildNiceOPdict:
    @pytest.fixture
    def systems(self):
        from fairmd.lipids.core import initialize_databank

        return initialize_databank()

    def test_work_empty(self):
        from fairmd.lipids.auxiliary.opconvertor import build_nice_OPdict
        from fairmd.lipids.molecules import Lipid

        lipid = Lipid("POPE")
        lipid.register_mapping()
        rdict = build_nice_OPdict({}, lipid)
        assert rdict == {}

    def test_build_fragmented(self, systems):
        from fairmd.lipids.auxiliary.opconvertor import build_nice_OPdict
        from fairmd.lipids.api import get_OP

        sys = systems.loc(281)
        opdata = get_OP(sys)
        rdict = build_nice_OPdict(opdata["POPC"], sys.lipids["POPC"])
        assert isinstance(rdict, dict)  # dict expected
        assert "sn-1" in rdict
        assert "sn-2" in rdict  # fragments at the top level

    def test_cnames_pl(self, systems):
        from fairmd.lipids.auxiliary.opconvertor import build_nice_OPdict
        from fairmd.lipids.api import get_OP

        sys = systems.loc(281)
        opdata = get_OP(sys)
        rdict = build_nice_OPdict(opdata["POPC"], sys.lipids["POPC"])

        # C check numbers
        def has_c(cname: str, flist: dict) -> bool:
            return any(_c["C"] == cname for _c in flist)

        check.is_true(has_c("2", rdict["sn-1"]))
        check.is_true(has_c("16", rdict["sn-1"]))
        check.is_false(has_c("1", rdict["sn-1"]))
        check.is_false(has_c("17", rdict["sn-1"]))
        # H check names
        check.is_true(all(_c["H"] in ["1", "2", "3"] for _c in rdict["sn-1"]))
        check.is_true(all(_c["H"] in ["1", "2", "3"] for _c in rdict["sn-2"]))
        # check backbone
        check.is_true(has_c("g1", rdict["glycerol backbone"]))
        check.is_true(has_c("g2", rdict["glycerol backbone"]))
        check.is_true(has_c("g3", rdict["glycerol backbone"]))
