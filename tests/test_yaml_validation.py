import copy
import os

import pytest

pytestmark = [pytest.mark.adddata, pytest.mark.min]

valid = {
    "DOI": "10.5281/zenodo.11614468",
    "TRJ": "566.trj",
    "TPR": "566.tpr",
    "SOFTWARE": "gromacs",
    "PREEQTIME": 0,
    "TIMELEFTOUT": 0,
    "SYSTEM": "120POPC_8CHOL_3968SOL_303K",
    "SOFTWARE_VERSION": "5.0.4",
    "COMPOSITION": {
        "DOPC": {"NAME": "DOPC", "MAPPING": "mappingDOPCcharmm.yaml"},
        "SOL": {"NAME": "TIP3", "MAPPING": "mappingTIP3PCHARMMgui.yaml"},
    },
}


@pytest.fixture
def valid_instance():
    return copy.deepcopy(valid)


def test_valid(valid_instance):
    from fairmd.lipids.SchemaValidation.ValidateYAML import validate_info_dict

    errors = validate_info_dict(valid_instance)
    assert len(errors) == 0


def test_missing_required(valid_instance):
    from fairmd.lipids.SchemaValidation.ValidateYAML import validate_info_dict

    del valid_instance["DOI"]
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "required"


def test_wrong_type(valid_instance):
    from fairmd.lipids.SchemaValidation.ValidateYAML import validate_info_dict

    valid_instance["TRJ"] = 1
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "type"


def test_composition_extra_key(valid_instance):
    from fairmd.lipids.SchemaValidation.ValidateYAML import validate_info_dict

    valid_instance["COMPOSITION"]["DOPC"]["NONSENSE"] = 1
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "additionalProperties"


def test_composition_missing_mapping(valid_instance):
    from fairmd.lipids.SchemaValidation.ValidateYAML import validate_info_dict

    del valid_instance["COMPOSITION"]["DOPC"]["MAPPING"]
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "required"


def test_good_united_atom_dict(valid_instance):
    from fairmd.lipids.SchemaValidation.ValidateYAML import validate_info_dict

    valid_instance["UNITEDATOM_DICT"] = {"atom1": "oxygen", "atom2": "hydrogen"}
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 0


def test_united_atom_dict_wrong_type(valid_instance):
    from fairmd.lipids.SchemaValidation.ValidateYAML import validate_info_dict

    valid_instance["UNITEDATOM_DICT"] = {"atom1": "oxygen", "atomo2": 2}
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "type"


def test_united_atom_dict_is_null(valid_instance):
    from fairmd.lipids.SchemaValidation.ValidateYAML import validate_info_dict

    valid_instance["UNITEDATOM_DICT"] = None
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 0


def test_missing_tpr_non_gromacs(valid_instance):
    from fairmd.lipids.SchemaValidation.ValidateYAML import validate_info_dict

    valid_instance["SOFTWARE"] = "openMM"
    del valid_instance["TPR"]
    valid_instance["PDB"] = "pdb-test"
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 0


def test_missing_tpr_gromacs(valid_instance):
    from fairmd.lipids.SchemaValidation.ValidateYAML import validate_info_dict

    del valid_instance["TPR"]
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "required"


def test_valid_info_file():
    from fairmd.lipids import FMDL_DATA_PATH
    from fairmd.lipids.SchemaValidation.ValidateYAML import validate_info_file

    valid_info_path = os.path.join(FMDL_DATA_PATH, "info", "info566.yaml")
    errors = validate_info_file(valid_info_path)
    assert len(errors) == 0
