"""Test the Molecules class for correct lipid metadata."""

import os
import pytest
import pytest_check as check

# run only on sim2 mocking data
pytestmark = [pytest.mark.rdkit]


def test_lipid_smarts():
    """Test atoms_by() method with SMARTS queries."""
    from fairmd.lipids.molecules import Lipid, MoleculeError

    popc = Lipid("POPC")
    _ = popc.rdkit_object  # rdkit object can be created without mapping

    _save = popc._metadata["bioschema_properties"]["smiles"]
    del popc._metadata["bioschema_properties"]["smiles"]
    with check.raises(MoleculeError):
        _ = popc.rdkit_object  # rdkit object cannot be created without smiles
    popc._metadata["bioschema_properties"]["smiles"] = _save

    with check.raises(ValueError):
        popc.atoms_by("~!x[", 0)

    with check.raises(KeyError):
        popc.atoms_by("N", -1)

    with check.raises(KeyError):
        popc.atoms_by("N", 1)  # 0 is max here

    popc.register_mapping("mappingPOPCcharmm.yaml")

    # Test that we can retrieve atoms by SMARTS
    nitrg_atoms = popc.atoms_by("N", 0)  # can call
    check.is_instance(nitrg_atoms, list, "Should return a list of atom universal names")
    check.equal(len(nitrg_atoms), 1, "Should have at least one nitrogen atom")
    check.equal(nitrg_atoms[0], "M_G3N6_M", "Bad atom selected by N from POPC")
