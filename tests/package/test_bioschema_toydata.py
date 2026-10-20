"""
System.bioschema over the ToyData simulations and their inherited template.

-------------------------------------------------------------------------------
NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import copy
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


def test_every_toy_system_stores_a_compact_block(systems):
    for system in systems:
        stored = system.readme["bioschema_properties"]
        assert not [k for k in stored["keywords"] if isinstance(k, dict)], system
        assert all(isinstance(t, str) for t in stored["measurementTechnique"]), system


def test_every_toy_system_resolves_the_inherited_terms(systems):
    validator = Draft7Validator(block_schema())
    for system in systems:
        stored = copy.deepcopy(system.readme["bioschema_properties"])
        block = system.bioschema

        codes = [k["termCode"] for k in block["keywords"] if isinstance(k, dict)]
        assert codes == INHERITED_KEYWORDS, system
        assert block["keywords"][len(codes) :] == stored["keywords"], system
        assert block["measurementTechnique"][0]["termCode"] == "operation_2476", system
        assert block["name"] == stored["name"], system
        assert not list(validator.iter_errors(block)), system

        # Resolving never writes into the README a writer would save back.
        assert system.readme["bioschema_properties"] == stored, system


def test_a_system_without_a_block_has_no_resolved_one(systems):
    from fairmd.lipids.core import System

    system = next(iter(systems))
    readme = {key: value for key, value in system.readme.items() if key != "bioschema_properties"}
    assert System(readme).bioschema is None
