"""
System.bioschema over the ToyData simulations and their inherited template.

-------------------------------------------------------------------------------
NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import copy
import os
import importlib.resources
import json

import pytest
from jsonschema import Draft7Validator

pytestmark = [pytest.mark.sim2, pytest.mark.min]

INHERITED_KEYWORDS = ["topic_3892", "topic_0176", "topic_0820", "topic_3306", "data_3870"]


def block_schema():
    schema_dir = importlib.resources.files("fairmd.lipids.schema_validation.schema")
    schema = json.loads((schema_dir / "simulation_schema.json").read_text(encoding="utf-8"))
    block = schema["properties"]["bioschema_properties"]
    block["definitions"] = schema["definitions"]
    return block


@pytest.fixture(scope="module")
def systems():
    from fairmd.lipids.core import initialize_databank

    return initialize_databank()


def test_a_system_is_a_metadata_provider(systems):
    from fairmd.lipids import FMDL_SIMU_PATH
    from fairmd.lipids._base import MetadataProvider
    from fairmd.lipids.core import System

    assert System.bioschema_template_dir() == FMDL_SIMU_PATH
    for system in systems:
        assert isinstance(system, MetadataProvider)


def test_every_public_view_of_the_block_is_the_resolved_one(systems):
    """No accessor outside the class returns the compact block the README stores."""
    for system in systems:
        block = system.bioschema
        assert system.metadata["bioschema_properties"] == block, system
        assert system["bioschema_properties"] == block, system
        assert system.get("bioschema_properties") == block, system
        assert system.readme["bioschema_properties"] == block, system
        assert dict(system)["bioschema_properties"] == block, system
        assert dict(system.items())["bioschema_properties"] == block, system


def test_every_toy_system_stores_a_compact_block(systems):
    for system in systems:
        for stored in (system._store["bioschema_properties"], system.to_readme()["bioschema_properties"]):
            assert not [k for k in stored["keywords"] if isinstance(k, dict)], system
            assert all(isinstance(t, str) for t in stored["measurementTechnique"]), system


def test_every_toy_system_resolves_the_inherited_terms(systems):
    validator = Draft7Validator(block_schema())
    for system in systems:
        stored = copy.deepcopy(system._store["bioschema_properties"])
        block = system.bioschema

        codes = [k["termCode"] for k in block["keywords"] if isinstance(k, dict)]
        assert codes == INHERITED_KEYWORDS, system
        assert block["keywords"][len(codes) :] == stored["keywords"], system
        assert block["measurementTechnique"][0]["termCode"] == "operation_2476", system
        assert block["name"] == stored["name"], system
        assert not list(validator.iter_errors(block)), system

        # Resolving never writes into what the README stores.
        assert system._store["bioschema_properties"] == stored, system


def test_writing_a_readme_back_keeps_the_block_compact(systems, tmp_path):
    """What match_experiments does: to_readme() is the stored form, not the resolved one."""
    import yaml

    from fairmd.lipids import FMDL_SIMU_PATH
    from fairmd.lipids.auxiliary.yaml_format import encode_canonical_yaml

    for system in systems:
        on_disk = yaml.safe_load(open(os.path.join(FMDL_SIMU_PATH, system["path"], "README.yaml"), encoding="utf-8"))[
            "bioschema_properties"
        ]
        out = tmp_path / "README.yaml"
        out.write_text(encode_canonical_yaml(system.to_readme()), encoding="utf-8")
        assert yaml.safe_load(out.read_text(encoding="utf-8"))["bioschema_properties"] == on_disk, system

        # A copy: changing it does not reach the record.
        system.to_readme()["bioschema_properties"]["keywords"].append("leaked")
        assert "leaked" not in system._store["bioschema_properties"]["keywords"], system


def test_assigning_the_resolved_block_back_stores_it_compact(systems):
    for system in copy.deepcopy(list(systems)):
        stored = copy.deepcopy(system._store["bioschema_properties"])
        system["bioschema_properties"] = system["bioschema_properties"]
        assert system._store["bioschema_properties"] == stored, system


def test_the_public_views_are_read_only(systems):
    """A change to a view would be lost, so it raises instead."""
    system = next(s for s in systems if "license" in s["bioschema_properties"])
    view = system["bioschema_properties"]
    for change in (
        lambda: view.__setitem__("name", "new"),
        lambda: view["keywords"].append("X"),
        lambda: view["license"].update(spdx="CC0-1.0"),
        lambda: system.metadata.__setitem__("TEMPERATURE", -1),
        lambda: system.metadata["bioschema_properties"]["keywords"].append("X"),
        lambda: system.readme.pop("ID"),
    ):
        with pytest.raises(TypeError, match="read-only"):
            change()
    # Values other than the block are the record's own, as system[key] gives them.
    assert system.metadata["COMPOSITION"] is system._store["COMPOSITION"]


def test_bioschema_and_copies_are_plain_mutable_values(systems):
    system = next(iter(systems))
    for value in (system.bioschema, copy.deepcopy(system["bioschema_properties"])):
        assert type(value) is dict and type(value["keywords"]) is list
        value["keywords"].append("mine")
    assert "mine" not in system["bioschema_properties"]["keywords"]


def test_a_view_serialises_to_json_but_cannot_be_dumped_as_a_readme(systems):
    import json

    import yaml

    system = next(iter(systems))
    assert json.loads(json.dumps(system["bioschema_properties"])) == system.bioschema
    # The guard against baking the template into a file: writers use to_readme().
    with pytest.raises(yaml.representer.RepresenterError):
        yaml.safe_dump(system.metadata)
    yaml.safe_dump(system.to_readme())


def test_a_system_built_from_a_resolved_view_stores_the_block_compact(systems):
    from fairmd.lipids.core import System

    for system in systems:
        rebuilt = System(dict(system))
        assert rebuilt.to_readme()["bioschema_properties"] == system.to_readme()["bioschema_properties"], system
        assert rebuilt.bioschema == system.bioschema, system


def test_the_resolved_block_is_cached_until_the_stored_one_changes(systems):
    system = copy.deepcopy(next(iter(systems)))
    first = system["bioschema_properties"]
    assert system["bioschema_properties"] is first
    assert system.metadata["bioschema_properties"] is first

    block = system.bioschema
    block["name"] = "renamed"
    system["bioschema_properties"] = block
    assert system["bioschema_properties"] is not first
    assert system["bioschema_properties"]["name"] == "renamed"


def test_a_system_without_a_block_has_no_resolved_one(systems):
    from fairmd.lipids.core import System

    system = next(iter(systems))
    readme = {key: value for key, value in system.to_readme().items() if key != "bioschema_properties"}
    bare = System(readme)
    assert bare.bioschema is None
    assert "bioschema_properties" not in bare.metadata
    assert "bioschema_properties" not in bare
