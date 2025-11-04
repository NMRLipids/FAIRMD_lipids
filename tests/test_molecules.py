"""Test the Molecules class for correct lipid metadata."""

import pytest
import pytest_check as check

# run only on sim2 mocking data
pytestmark = [pytest.mark.sim2, pytest.mark.min]

LIPIDS_SET_LENGTH = 5
POPE_MOLECULAR_WEIGHT = 718


def test_lipids_metadata():
    """Test metadata of lipids_set, especially for POPE."""
    from fairmd.lipids import FMDL_MOL_PATH
    from fairmd.lipids.settings.molecules import lipids_set

    print(FMDL_MOL_PATH)
    print(lipids_set)

    check.equal(len(lipids_set), LIPIDS_SET_LENGTH, "LipidSet should have length 5")
    check.is_in("POPE", lipids_set, "POPE should be in lipids_set")
    pope = lipids_set.get("POPE")
    assert pope is not None, "POPE should not be None"
    check.equal(pope.name, "POPE", "POPE name should be 'POPE'")
    print(f"POPE metadata: {pope.metadata}")
    metadata = pope.metadata
    assert metadata is not None, "POPE metadata should not be None"
    check.is_not_none(metadata.get("bioschema_properties"), "POPE bioschema_properties should not be None")
    check.equal(metadata["NMRlipids"]["id"], "POPE", "POPE bioschema_properties name should be 'POPE'")
    check.equal(
        metadata["NMRlipids"]["name"],
        "1-palmitoyl-2-oleoyl-sn-glycero-3-phosphoethanolamine",
        "POPE NMRlipids name should match expected value",
    )
    check.equal(
        metadata["bioschema_properties"]["molecularWeight"],
        POPE_MOLECULAR_WEIGHT,
        "POPE bioschema_properties molecularWeight should be 718.0",
    )
    check.equal(
        metadata["bioschema_properties"]["molecularFormula"],
        "C39H76NO8P",
        "POPE bioschema_properties formula should be 'C39H76NO8P'",
    )
    check.equal(
        metadata["bioschema_properties"]["smiles"],
        "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OCCN)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
        "POPE bioschema_properties smiles should match expected value",
    )
