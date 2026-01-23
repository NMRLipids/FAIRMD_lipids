"""Test the Molecules class for correct lipid metadata."""

import os
import pytest
import pytest_check as check

# run only on sim2 mocking data
pytestmark = [pytest.mark.sim2, pytest.mark.min]

LIPIDS_SET_LENGTH = 5
POPE_MOLECULAR_WEIGHT = 717.5


@pytest.fixture(scope="module")
def toy_mols() -> dict:
    from fairmd.lipids import FMDL_MOL_PATH
    from fairmd.lipids.molecules import lipids_set

    print(FMDL_MOL_PATH)
    print(lipids_set)

    mol1 = lipids_set.get("POPE")
    mol1.register_mapping("mappingPOPEcharmm.yaml")

    mol2 = lipids_set.get("POPC")
    mol2.register_mapping("mappingPOPClipid14.yaml")

    return {"pope/charmm": mol1, "popc/amber": mol2}


def test_mapping_dict():
    from fairmd.lipids.molecules import lipids_set, MoleculeError

    mol1 = lipids_set.get("POPE")
    with pytest.raises(MoleculeError, match="not registered"):
        _ = mol1.mapping_dict


def test_uan2selection(toy_mols):
    toy_pope = toy_mols["pope/charmm"]
    selstr = toy_pope.uan2selection("M_G1C2_M", "POPE")
    check.equal(selstr, "name C31 and resname POPE", "Selection string should match expected value (charmm36)")
    with check.raises(KeyError):
        toy_pope.uan2selection("NON_EXISTENT_UAN", "POPE")

    toy_popc = toy_mols["popc/amber"]
    selstr = toy_popc.uan2selection("M_G1C6H1_M", "POPC")
    check.equal(selstr, "name H5R and resname PA", "Selection string should match expected value (lipid14)")
    with check.raises(KeyError):
        toy_pope.uan2selection("NON_EXISTENT_UAN", "POPC")


def test_md2uan(toy_mols):
    import MDAnalysis as mda
    from fairmd.lipids.molecules import NonLipid, MoleculeMappingError

    popc2_fp = os.path.join(os.path.dirname(__file__), "misc_data", "popc2.gro")
    u_popc2 = mda.Universe(popc2_fp)
    toy_popc = toy_mols["popc/amber"]
    uan = toy_popc.md2uan("H5R", "PA")
    check.equal(uan, "M_G1C6H1_M", "Universal Atom Name should match expected value (lipid14)")
    uan = toy_popc.md2uan("H5R", "OL")
    check.equal(uan, "M_G2C6H1_M", "Universal Atom Name should match expected value (lipid14)")
    with check.raises(MoleculeMappingError):
        toy_popc.md2uan("NON_EXISTENT_MD_NAME", "PA")
    with check.raises(MoleculeMappingError):
        toy_popc.md2uan("H5R", "NON_EXISTENT_RESNAME")

    cl_mol = NonLipid("CLA")
    cl_mol.register_mapping("mappingCL.yaml")
    check.equal(cl_mol.md2uan("CL"), "M_Cl_M")
    check.equal(cl_mol.md2uan("Cl"), "M_Cl_M")
    with check.raises(MoleculeMappingError):
        cl_mol.md2uan("Cl-")


def test_check_mapping_amber(toy_mols):
    import MDAnalysis as mda
    from fairmd.lipids.molecules import MoleculeMappingError

    popc2_fp = os.path.join(os.path.dirname(__file__), "misc_data", "popc2.gro")
    popc1_fp = os.path.join(os.path.dirname(__file__), "misc_data", "popc1.gro")

    u_popc2 = mda.Universe(popc2_fp)  # correct amber

    toy_popc = toy_mols["popc/amber"]

    # correct mapping
    with check.check:
        toy_popc.check_mapping(u_popc2, "POPC")

    # correct mapping, incorrect resname
    # it doesn't matter here because resnames are in the mapping file
    with check.check:
        toy_popc.check_mapping(u_popc2, "BUBA")

    # mapping is OK but one atom in the mol is missing
    u_popc2_minus_atom = mda.Merge(u_popc2.select_atoms("id 2:1000"))
    with check.raises(MoleculeMappingError):
        toy_popc.check_mapping(u_popc2_minus_atom, "")

    # incorrect residue name in MD
    u_popc2.select_atoms("resname PA").residues.resnames = "PAA"
    with check.raises(MoleculeMappingError):
        toy_popc.check_mapping(u_popc2, "")


def test_check_mapping(toy_mols):
    import MDAnalysis as mda
    from fairmd.lipids.molecules import MoleculeMappingError

    pope1_fp = os.path.join(os.path.dirname(__file__), "misc_data", "pope1.gro")
    popc1_fp = os.path.join(os.path.dirname(__file__), "misc_data", "popc1.gro")

    u_pope = mda.Universe(pope1_fp)
    u_popc = mda.Universe(popc1_fp)

    toy_pope = toy_mols["pope/charmm"]

    # correct mapping
    with check.check:
        toy_pope.check_mapping(u_pope, "POPE")

    # correct mapping, incorrect resname
    with check.raises(MoleculeMappingError):
        toy_pope.check_mapping(u_pope, "BUBA")

    # correct resname, incorrect mapping
    with check.raises(MoleculeMappingError):
        toy_pope.check_mapping(u_popc, "POPC")

    # mapping is OK but one atom in the mol is missing
    u_pope_minus_atom = mda.Merge(u_pope.select_atoms("id 2:1000"))
    with check.raises(MoleculeMappingError):
        toy_pope.check_mapping(u_pope_minus_atom, "POPE")


def test_lipids_metadata():
    """Test metadata of lipids_set, especially for POPE."""
    from fairmd.lipids import FMDL_MOL_PATH
    from fairmd.lipids.molecules import lipids_set

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
        "POPE bioschema_properties molecularWeight should be 718.5",
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
