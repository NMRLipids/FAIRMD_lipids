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


def only_error(errors):
    """The single validation error, named by the property it is about.

    Returns ``(dotted path, validator)``. Asserting on the validator alone lets
    an unrelated failure satisfy a test -- a block can be rejected by `required`
    for a reason that has nothing to do with the field the test removed.
    """
    assert len(errors) == 1, [f"{list(e.absolute_path)}: {e.message}" for e in errors]
    error = errors[0]
    return ".".join(str(part) for part in error.absolute_path), error.validator


bioschema = {
    # As the generator writes it: composed from the entry's own fields, never the
    # deposition title.
    "name": (
        "Molecular dynamics trajectory of a POPC bilayer at 313 K "
        "(CHARMM36, GROMACS 5.0.4, 200 ns)"
    ),
    "alternateName": "200POPC_9000SOL_313K",
    "description": (
        "Molecular dynamics trajectory of a lipid bilayer of 200 POPC at 313 K. "
        "Simulated with CHARMM36 in GROMACS 5.0.4 for 200 ns (40000 atoms). "
        "Part of the NMRlipids Databank and available from "
        "https://doi.org/10.5281/zenodo.4040423."
    ),
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
    "isPartOf": {
        "@type": "Dataset",
        "@id": "https://doi.org/10.5281/zenodo.4040423",
        "identifier": "10.5281/zenodo.4040423",
        "url": "https://doi.org/10.5281/zenodo.4040423",
        "name": "Simulation trajectories of POPC bilayers",
        "publisher": "Zenodo",
    },
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
    assert only_error(errors) == ("bioschema_properties", "additionalProperties")
    assert "nonsense" in errors[0].message


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
    assert only_error(errors) == ("bioschema_properties", "required")
    assert "'name'" in errors[0].message


def test_readme_bioschema_creator_missing_name(readme_with_bioschema):
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    readme_with_bioschema["bioschema_properties"]["creator"].append({"identifier": "x"})
    errors = validate_readme_dict(readme_with_bioschema)
    assert only_error(errors) == ("bioschema_properties.creator.1", "required")
    assert "'name'" in errors[0].message


def test_readme_bioschema_bad_date(readme_with_bioschema):
    """datePublished is YYYY, YYYY-MM or YYYY-MM-DD."""
    from fairmd.lipids.schema_validation.validate_yaml import validate_readme_dict

    readme_with_bioschema["bioschema_properties"]["datePublished"] = "21/09/2020"
    errors = validate_readme_dict(readme_with_bioschema)
    assert only_error(errors) == ("bioschema_properties.datePublished", "pattern")


# ---------------------------------------------------------------------------
# experiment_schema.json
#
# Experiment READMEs have no validate_* entry point yet, so the schema is driven
# directly. The DOI rules are the part contributors get wrong, and the part the
# documentation and the enrichment tooling both depend on, so they are pinned
# here rather than left to the generator tests.
# ---------------------------------------------------------------------------

ARTICLE_DOI = "10.1039/B508190D"
DATA_DOI = "10.18710/ETWNCU"


def experiment_schema():
    import json

    from fairmd.lipids.schema_validation import validate_yaml

    path = os.path.join(os.path.dirname(validate_yaml.__file__), "schema", "experiment_schema.json")
    with open(path) as handle:
        return json.load(handle)


def experiment_errors(instance):
    from jsonschema import Draft7Validator

    return sorted(Draft7Validator(experiment_schema()).iter_errors(instance), key=lambda e: e.path)


@pytest.fixture
def valid_experiment():
    return {
        "ARTICLE_DOI": ARTICLE_DOI,
        "TEMPERATURE": 314.0,
        "MEMBRANE_COMPOSITION": {"DMPC": 1},
        "TOTAL_HYDRATION": 95,
        "PH": "UNKNOWN",
        "NMR": {
            "METHOD": "2H:QE",
            "INSTRUMENT": "Bruker MSL-300",
            "SIGN_MEASURED": "NONE",
            "T_RF_HEATING": "UNKNOWN",
        },
    }


def test_valid_experiment(valid_experiment):
    assert experiment_errors(valid_experiment) == []


def test_experiment_with_both_dois_is_valid(valid_experiment):
    """Both may be given; which one gets cited is the generator's business."""
    valid_experiment["DATA_DOI"] = DATA_DOI
    assert experiment_errors(valid_experiment) == []


def test_experiment_with_only_a_data_doi_is_valid(valid_experiment):
    del valid_experiment["ARTICLE_DOI"]
    valid_experiment["DATA_DOI"] = DATA_DOI
    assert experiment_errors(valid_experiment) == []


def test_experiment_without_any_doi_is_rejected(valid_experiment):
    del valid_experiment["ARTICLE_DOI"]
    errors = experiment_errors(valid_experiment)
    assert len(errors) == 1
    assert errors[0].validator == "anyOf"
    assert "ARTICLE_DOI" in str(errors[0]) and "DATA_DOI" in str(errors[0])


def test_data_ref_does_not_stand_in_for_a_doi(valid_experiment):
    """DATA_REF is free text, so it cannot satisfy the requirement on its own."""
    del valid_experiment["ARTICLE_DOI"]
    valid_experiment["DATA_REF"] = "nmrXiv sample S-1234, no DOI assigned"
    errors = experiment_errors(valid_experiment)
    assert len(errors) == 1
    assert errors[0].validator == "anyOf"


@pytest.mark.parametrize(
    "value",
    [
        "https://doi.org/10.1039/B508190D",
        "doi:10.1039/B508190D",
        "10.1039",
        "10.1039/",
        "",
    ],
)
def test_dois_must_be_bare(valid_experiment, value):
    """Both fields take 10.xxxx/yyyy: no URL, no prefix, nothing empty."""
    valid_experiment["ARTICLE_DOI"] = value
    paths = {".".join(str(p) for p in e.absolute_path) for e in experiment_errors(valid_experiment)}
    assert paths == {"ARTICLE_DOI"}


def test_experiment_template_covers_every_schema_field():
    """The template a contributor copies must not drift from what is accepted."""
    import yaml as yaml_module

    from fairmd.lipids.schema_validation import validate_yaml

    path = os.path.join(os.path.dirname(validate_yaml.__file__), "schema", "experiment_template.yaml")
    with open(path) as handle:
        template = yaml_module.safe_load(handle)

    schema = experiment_schema()
    # bioschema_properties is written by the enrichment tooling, so it is
    # deliberately absent from a template meant to be filled in by hand.
    assert set(schema["properties"]) - set(template) == {"bioschema_properties"}
    assert set(template) - set(schema["properties"]) == set()
    for block in ("NMR", "XRAY"):
        assert set(template[block]) == set(schema["properties"][block]["properties"])
