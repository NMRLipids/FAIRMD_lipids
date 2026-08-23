"""
Unit-testing for NMRPCA analysis module

Test data is stored in `./ToyData/Simulations.1`

NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import glob
import logging
import os

import pytest
import pytest_check as check

pytestmark = pytest.mark.sim1


def test_parser_verify_lipid():
    import fairmd.lipids as fmdl
    from fairmd.lipids.analib.analyze_nmrpca import Parser
    from fairmd.lipids.core import lipids_set, Lipid

    # Load databank to get lipid definitions
    popc: Lipid = lipids_set.get("POPC")
    popc.register_mapping("mappingPOPCcharmm.yaml")
    check.is_true(Parser.verify_lipid(popc))
    popc.register_mapping("mappingPOPClipid14.yaml")
    check.is_true(Parser.verify_lipid(popc))

    cholesterol: Lipid = lipids_set.get("CHOL")
    cholesterol.register_mapping("mappingCHOLESTEROLlipid14.yaml")
    check.is_false(Parser.verify_lipid(cholesterol))

    tocl: Lipid = lipids_set.get("TOCL")
    tocl.register_mapping("mappingTOCLcharmm.yaml")
    check.is_false(Parser.verify_lipid(tocl))

    pope: Lipid = lipids_set.get("POPE")
    pope.register_mapping("mappingPOPEGROMOS43A1-S3.yaml")
    check.is_true(Parser.verify_lipid(pope))
