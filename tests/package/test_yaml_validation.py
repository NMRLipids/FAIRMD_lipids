import copy
import os

import pytest

pytestmark = [pytest.mark.adddata, pytest.mark.min]

valid = {
    "DOI": "10.5281/zenodo.11614468",
    "TRJ": "566.trj",
    "TPR": "566.tpr",
    "SOFTWARE": "gromacs",
    "PREEQTIME": 2.0,
    "TIMELEFTOUT": 10,
    "SYSTEM": "120POPC_8CHOL_3968SOL_303K",
    "SOFTWARE_VERSION": "5.0.4",
    "FF": "CHARMM36",
    "AUTHORS_CONTACT": "Einstein, Albert",
    "COMPOSITION": {
        "DOPC": {"NAME": "DOPC", "MAPPING": "mappingDOPCcharmm.yaml"},
        "SOL": {"NAME": "TIP3", "MAPPING": "mappingTIP3PCHARMMgui.yaml"},
    },
}


@pytest.fixture(scope="module")
def systems():
    from fairmd.lipids import FMDL_DATA_PATH, FMDL_SIMU_PATH


@pytest.fixture
def valid_instance():
    return copy.deepcopy(valid)


def test_valid(valid_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_info_dict

    errors = validate_info_dict(valid_instance)
    assert len(errors) == 0


def test_missing_required(valid_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_info_dict

    del valid_instance["DOI"]
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "required"


def test_wrong_type(valid_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_info_dict

    valid_instance["TRJ"] = 1
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "type"


def test_composition_extra_key(valid_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_info_dict

    valid_instance["COMPOSITION"]["DOPC"]["NONSENSE"] = 1
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "additionalProperties"


def test_composition_missing_mapping(valid_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_info_dict

    del valid_instance["COMPOSITION"]["DOPC"]["MAPPING"]
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "required"


def test_good_united_atom_dict(valid_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_info_dict

    valid_instance["UNITEDATOM_DICT"] = {"atom1": "oxygen", "atom2": "hydrogen"}
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 0


def test_united_atom_dict_wrong_type(valid_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_info_dict

    valid_instance["UNITEDATOM_DICT"] = {"atom1": "oxygen", "atomo2": 2}
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "type"


def test_united_atom_dict_is_null(valid_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_info_dict

    valid_instance["UNITEDATOM_DICT"] = None
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 0


def test_missing_tpr_non_gromacs(valid_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_info_dict

    valid_instance["SOFTWARE"] = "openMM"
    del valid_instance["TPR"]
    valid_instance["PDB"] = "test.pdb"
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 0


def test_missing_tpr_gromacs(valid_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_info_dict

    del valid_instance["TPR"]
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "required"


def test_missing_FF(valid_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_info_dict

    del valid_instance["FF"]
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "required"


def test_missing_authors_contact(valid_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_info_dict

    del valid_instance["AUTHORS_CONTACT"]
    errors = validate_info_dict(valid_instance)
    assert len(errors) == 1
    assert errors[0].validator == "required"


def test_valid_info_file():
    from fairmd.lipids import FMDL_DATA_PATH
    from fairmd.lipids.schema_validation.validate_yaml import validate_info_file

    valid_info_path = os.path.join(FMDL_DATA_PATH, "info", "info566.yaml")
    errors = validate_info_file(valid_info_path)
    assert len(errors) == 0


@pytest.fixture
def valid_readme_instance(valid_instance):
    """Base info dict + README-specific required fields and changes."""
    inst = copy.deepcopy(valid_instance)

    inst["TRJ"] = [[inst["TRJ"]]]
    inst["TPR"] = [[inst["TPR"]]]

    inst.update(
        {
            "TRAJECTORY_SIZE": 123456,
            "TRJLENGTH": 1000.0,
            "NUMBER_OF_ATOMS": 12345,
            "DATEOFRUNNING": "2024-01-01",
            "ID": 1,
        }
    )

    for comp in inst["COMPOSITION"].values():
        comp["COUNT"] = 1
    return inst


@pytest.fixture
def valid_readme_namd(valid_readme_instance):
    inst = copy.deepcopy(valid_readme_instance)
    inst["SOFTWARE"] = "NAMD"
    del inst["TPR"]
    inst["PSF"] = [["valid.psf"]]
    return inst


@pytest.fixture
def valid_readme_openMM(valid_readme_instance):
    inst = copy.deepcopy(valid_readme_instance)
    inst["SOFTWARE"] = "openMM"
    inst["AMBERTOP"] = [["valid.top"]]
    return inst


def test_valid_readme(valid_readme_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    errors = validate_readme_dict(valid_readme_instance)
    assert len(errors) == 0


def test_readme_missing_required(valid_readme_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    del valid_readme_instance["TRAJECTORY_SIZE"]
    errors = validate_readme_dict(valid_readme_instance)
    assert len(errors) == 1
    assert errors[0].validator == "required"


def test_readme_wrong_type_id(valid_readme_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    valid_readme_instance["ID"] = "not-an-int"
    errors = validate_readme_dict(valid_readme_instance)
    assert len(errors) == 1
    assert errors[0].validator == "type"


def test_readme_wrong_type_traj_size(valid_readme_instance):
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    valid_readme_instance["TRAJECTORY_SIZE"] = -1
    errors = validate_readme_dict(valid_readme_instance)
    assert len(errors) == 1
    assert errors[0].validator in ("minimum", "type")


def test_valid_readme_namd(valid_readme_namd):
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    errors = validate_readme_dict(valid_readme_namd)
    assert len(errors) == 0


def test_namd_wrong_file_ending(valid_readme_namd):
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    valid_readme_namd["PSF"] = ["bad.json"]

    errors = validate_readme_dict(valid_readme_namd)
    assert len(errors) == 1


def test_valid_readme_openNN(valid_readme_openMM):
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    errors = validate_readme_dict(valid_readme_openMM)
    assert len(errors) == 0


def test_wrong_filetype_readme_openNN(valid_readme_openMM):
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    valid_readme_openMM["AMBERTOP"] = "bad.gro"
    errors = validate_readme_dict(valid_readme_openMM)
    assert len(errors) == 1


def test_valid_readme_file():
    from fairmd.lipids import FMDL_DATA_PATH
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_file

    valid_info_path = os.path.join(FMDL_DATA_PATH, "Simulations.2", "aa0", "README.yaml")
    errors = validate_readme_file(valid_info_path)
    assert len(errors) == 0


bioschema = {
    "name": "POPC CHARMM36 T313K",
    "alternateName": "200POPC_9000SOL_313K",
    "description": "A POPC bilayer of 200 lipids simulated at 313 K.",
    "sameAs": "https://doi.org/10.5281/zenodo.4040423",
    "datePublished": "2020-09-21",
    "license": {
        "spdx": "CC-BY-4.0",
        "name": "Creative Commons Attribution 4.0 International",
        "url": "https://spdx.org/licenses/CC-BY-4.0.html",
        "sameAs": "https://creativecommons.org/licenses/by/4.0/legalcode",
    },
    "publisher": "Zenodo",
    "creator": [{"name": "Ollila", "identifier": "https://orcid.org/0000-0002-8135-3562"}],
    "citation": ["10.1021/acs.jctc.5b00935"],
    "keywords": [
        {
            "@type": "DefinedTerm",
            "name": "Molecular dynamics",
            "termCode": "topic_0176",
            "inDefinedTermSet": "http://edamontology.org",
            "url": "http://edamontology.org/topic_0176",
        },
        "POPC",
    ],
    "measurementTechnique": ["Molecular dynamics simulation (gromacs 5)"],
    "variableMeasured": ["area per lipid"],
    "distribution": [
        {
            "@type": "DataDownload",
            "contentUrl": "https://doi.org/10.5281/zenodo.4040423",
            "name": "200POPC_9000SOL_313K",
            "contentSize": 3979765856,
            "encodingFormat": "application/x-xtc",
            "hasPart": ["run.xtc", "run.tpr"],
        }
    ],
    "isBasedOn": ["experiments/FormFactors/10.1016/j.bbamem.2011.07.022/11"],
    "accessRights": "openAccess",
    "_source": {
        "api": "datacite",
        "doi": "10.5281/zenodo.4040423",
        "spdxLicenseList": "3.28.0",
        "retrieved": "2026-09-08",
    },
}


@pytest.fixture
def readme_with_bioschema(valid_readme_instance):
    inst = copy.deepcopy(valid_readme_instance)
    inst["bioschema_properties"] = copy.deepcopy(bioschema)
    return inst


def test_readme_with_bioschema_properties(readme_with_bioschema):
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    errors = validate_readme_dict(readme_with_bioschema)
    assert len(errors) == 0


def test_readme_without_bioschema_properties(valid_readme_instance):
    """The block is optional: entries that are not enriched yet stay valid."""
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    errors = validate_readme_dict(valid_readme_instance)
    assert len(errors) == 0


def test_readme_bioschema_unknown_property(readme_with_bioschema):
    """Only Bioschemas Dataset profile properties are accepted."""
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    readme_with_bioschema["bioschema_properties"]["nonsense"] = 1
    errors = validate_readme_dict(readme_with_bioschema)
    assert len(errors) == 1
    assert errors[0].validator == "additionalProperties"


def test_readme_bioschema_dct_prefixed_allowed(readme_with_bioschema):
    """DCMI terms are pre-approved, as the profile itself does for dct:conformsTo."""
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    readme_with_bioschema["bioschema_properties"]["dct:accessRights"] = "openAccess"
    errors = validate_readme_dict(readme_with_bioschema)
    assert len(errors) == 0


def test_readme_bioschema_missing_name(readme_with_bioschema):
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    del readme_with_bioschema["bioschema_properties"]["name"]
    errors = validate_readme_dict(readme_with_bioschema)
    assert len(errors) == 1
    assert errors[0].validator == "required"


def test_readme_bioschema_creator_missing_name(readme_with_bioschema):
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    readme_with_bioschema["bioschema_properties"]["creator"].append({"identifier": "x"})
    errors = validate_readme_dict(readme_with_bioschema)
    assert len(errors) == 1
    assert errors[0].validator == "required"


def test_readme_bioschema_bad_date(readme_with_bioschema):
    """datePublished is YYYY, YYYY-MM or YYYY-MM-DD."""
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    readme_with_bioschema["bioschema_properties"]["datePublished"] = "21/09/2020"
    errors = validate_readme_dict(readme_with_bioschema)
    assert len(errors) == 1
    assert errors[0].validator == "pattern"
