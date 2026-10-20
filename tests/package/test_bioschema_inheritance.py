"""
The inheritance of ``bioschema_properties`` from the templates above a record.

-------------------------------------------------------------------------------
NOTE: globally import of fairmd-lipids is **STRICTLY FORBIDDEN** because it
      breaks the substitution of global path folders
"""

import copy
import os

import pytest
import yaml

pytestmark = [pytest.mark.nodata, pytest.mark.min]

TEMPLATE_FILENAME = "bioschema_properties.yaml"


@pytest.fixture
def rules():
    from fairmd.lipids import bioschema_inheritance

    return bioschema_inheritance


MD_TOPIC = {"@type": "DefinedTerm", "name": "Molecular dynamics", "termCode": "topic_0176"}
BIOPHYSICS = {"@type": "DefinedTerm", "name": "Biophysics", "termCode": "topic_3306"}
CC_BY = {
    "spdx": "CC-BY-4.0",
    "name": "Creative Commons Attribution 4.0 International",
    "url": "https://spdx.org/licenses/CC-BY-4.0.html",
    "sameAs": "https://creativecommons.org/licenses/by/4.0/",
}
CC0 = {"spdx": "CC0-1.0", "url": "https://spdx.org/licenses/CC0-1.0.html"}


def write_template(folder, content):
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / TEMPLATE_FILENAME
    path.write_text(yaml.safe_dump(content, sort_keys=False), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# The merge: a recursive dict update
# ---------------------------------------------------------------------------


def test_nested_dicts_are_merged_key_by_key_at_every_depth(rules):
    parent = {"includedInDataCatalog": {"@type": "DataCatalog", "provider": {"name": "NMRlipids", "url": "a"}}}
    child = {"includedInDataCatalog": {"provider": {"url": "b"}, "version": "2"}}
    merged = rules.merge_blocks(parent, child)
    assert merged["includedInDataCatalog"] == {
        "@type": "DataCatalog",
        "provider": {"name": "NMRlipids", "url": "b"},
        "version": "2",
    }


def test_a_closer_scalar_wins(rules):
    assert rules.merge_blocks({"publisher": "Zenodo"}, {"publisher": "nmrXiv"}) == {"publisher": "nmrXiv"}


def test_a_type_mismatch_is_won_by_the_closer_value(rules):
    assert rules.merge_blocks({"version": {"major": 1}}, {"version": "1.0"}) == {"version": "1.0"}


def test_lists_are_extended_inherited_items_first_without_duplicates(rules):
    merged = rules.merge_blocks({"keywords": [MD_TOPIC, BIOPHYSICS]}, {"keywords": ["POPC", BIOPHYSICS, "POPS"]})
    assert merged["keywords"] == [MD_TOPIC, BIOPHYSICS, "POPC", "POPS"]


def test_an_inherited_property_the_record_does_not_state_is_added(rules):
    merged = rules.merge_blocks({"variableMeasured": ["order parameter"]}, {"name": "x"})
    assert merged == {"name": "x", "variableMeasured": ["order parameter"]}


def test_record_owned_properties_are_replaced_whole(rules):
    """Merged key by key, a CC0 record would keep the inherited CC-BY name and sameAs."""
    assert rules.REPLACED_WHOLE == {
        "license",
        "articleLicense",
        "isPartOf",
        "creator",
        "citation",
        "distribution",
        "isBasedOn",
    }
    for key in rules.REPLACED_WHOLE:
        parent_value = CC_BY if key in ("license", "articleLicense", "isPartOf") else ["inherited"]
        child_value = CC0 if isinstance(parent_value, dict) else ["own"]
        merged = rules.merge_blocks({key: parent_value}, {key: child_value})
        assert merged[key] == child_value, key


def test_the_licence_case_is_not_mixed(rules):
    merged = rules.merge_blocks({"license": CC_BY}, {"license": CC0})
    assert merged["license"] == CC0
    assert "name" not in merged["license"]
    assert "sameAs" not in merged["license"]


def test_record_identity_is_never_inherited(rules):
    assert rules.NOT_INHERITED == {"_source", "name", "alternateName", "description", "sameAs"}
    for key in rules.NOT_INHERITED:
        assert key not in rules.merge_blocks({key: "from a template"}, {}), key
        assert rules.merge_blocks({key: "from a template"}, {key: "own"})[key] == "own", key


def test_merging_leaves_both_inputs_untouched(rules):
    parent = {"keywords": [MD_TOPIC], "license": CC_BY, "catalog": {"a": {"b": 1}}}
    child = {"keywords": ["POPC"], "catalog": {"a": {"c": 2}}}
    parent_before, child_before = copy.deepcopy(parent), copy.deepcopy(child)

    merged = rules.merge_blocks(parent, child)
    merged["keywords"].append("changed")
    merged["license"]["spdx"] = "changed"
    merged["catalog"]["a"]["b"] = "changed"

    assert parent == parent_before
    assert child == child_before


# ---------------------------------------------------------------------------
# Compaction, the inverse used when a README is written
# ---------------------------------------------------------------------------

INHERITED = {
    "keywords": [MD_TOPIC, BIOPHYSICS],
    "license": CC_BY,
    "variableMeasured": [{"@type": "PropertyValue", "name": "C-H bond order parameter"}],
    "catalog": {"@type": "DataCatalog", "provider": {"name": "NMRlipids"}},
}

FULL = {
    "name": "C-H bond order parameters of a POPC bilayer",
    "keywords": [MD_TOPIC, BIOPHYSICS, "POPC"],
    "license": CC_BY,
    "variableMeasured": [{"@type": "PropertyValue", "name": "C-H bond order parameter"}],
    "catalog": {"@type": "DataCatalog", "provider": {"name": "NMRlipids", "url": "x"}},
    "publisher": "Elsevier",
}


def test_compaction_keeps_only_what_the_record_adds(rules):
    assert rules.compact_block(FULL, INHERITED) == {
        "name": "C-H bond order parameters of a POPC bilayer",
        "keywords": ["POPC"],
        "catalog": {"provider": {"url": "x"}},
        "publisher": "Elsevier",
    }


def test_a_record_value_that_differs_from_a_whole_one_is_kept_whole(rules):
    assert rules.compact_block({"license": CC0}, INHERITED) == {"license": CC0}


@pytest.mark.parametrize(
    "block",
    [
        FULL,
        {"name": "x", "license": CC0, "keywords": ["POPC", BIOPHYSICS]},
        {"name": "x", "catalog": "a string now", "variableMeasured": []},
        {"name": "x"},
    ],
)
def test_compaction_loses_nothing_the_resolved_block_states(rules, block):
    assert rules.merge_blocks(INHERITED, rules.compact_block(block, INHERITED)) == rules.merge_blocks(INHERITED, block)


def test_compacting_against_nothing_changes_nothing(rules):
    assert rules.compact_block(FULL, {}) == FULL


# ---------------------------------------------------------------------------
# The chain of folders
# ---------------------------------------------------------------------------


def test_templates_apply_outermost_first(rules, tmp_path):
    base = tmp_path / "experiments"
    write_template(base, {"license": CC_BY, "keywords": [BIOPHYSICS], "publisher": "outer"})
    write_template(base / "OrderParameters", {"keywords": [MD_TOPIC], "publisher": "inner"})
    record = base / "OrderParameters" / "10.1039" / "B508190D" / "1"
    record.mkdir(parents=True)

    assert rules.template_paths(record, base) == [
        (base / TEMPLATE_FILENAME).resolve(),
        (base / "OrderParameters" / TEMPLATE_FILENAME).resolve(),
    ]
    assert rules.inherited_block(record, base) == {
        "license": CC_BY,
        "keywords": [BIOPHYSICS, MD_TOPIC],
        "publisher": "inner",
    }
    resolved = rules.resolve_block({"name": "x", "keywords": ["DMPC"]}, record, base)
    assert resolved == {
        "name": "x",
        "keywords": [BIOPHYSICS, MD_TOPIC, "DMPC"],
        "license": CC_BY,
        "publisher": "inner",
    }


def test_a_template_in_the_record_folder_itself_is_not_inherited(rules, tmp_path):
    base = tmp_path / "Simulations"
    record = base / "aa0"
    write_template(record, {"publisher": "nobody"})
    assert rules.inherited_block(record, base) == {}


def test_a_template_above_the_base_is_not_inherited(rules, tmp_path):
    write_template(tmp_path, {"publisher": "nobody"})
    base = tmp_path / "Simulations"
    record = base / "aa0"
    record.mkdir(parents=True)
    assert rules.inherited_block(record, base) == {}


def test_a_record_outside_the_base_inherits_nothing(rules, tmp_path):
    base = tmp_path / "Simulations"
    write_template(base, {"publisher": "Zenodo"})
    elsewhere = tmp_path / "elsewhere" / "aa0"
    elsewhere.mkdir(parents=True)
    assert rules.inherited_block(elsewhere, base) == {}
    assert rules.resolve_block({"name": "x"}, elsewhere, base) == {"name": "x"}


def test_identity_keys_in_a_template_are_ignored(rules, tmp_path):
    base = tmp_path / "Simulations"
    write_template(base, {"name": "every simulation", "_source": {"api": "x"}, "publisher": "Zenodo"})
    record = base / "aa0"
    record.mkdir(parents=True)
    assert rules.inherited_block(record, base) == {"publisher": "Zenodo"}


def test_an_edited_template_is_read_again(rules, tmp_path):
    base = tmp_path / "Simulations"
    path = write_template(base, {"publisher": "first"})
    record = base / "aa0"
    record.mkdir(parents=True)
    assert rules.inherited_block(record, base) == {"publisher": "first"}

    path.write_text(yaml.safe_dump({"publisher": "second"}), encoding="utf-8")
    stat = path.stat()
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
    assert rules.inherited_block(record, base) == {"publisher": "second"}


def test_the_cached_template_cannot_be_changed_through_a_result(rules, tmp_path):
    base = tmp_path / "Simulations"
    path = write_template(base, {"keywords": [MD_TOPIC]})
    record = base / "aa0"
    record.mkdir(parents=True)

    rules.inherited_block(record, base)["keywords"].append("leaked")
    rules.load_template(path)["keywords"].append("leaked")
    assert rules.inherited_block(record, base) == {"keywords": [MD_TOPIC]}


def test_a_template_that_is_not_a_mapping_is_an_error(rules, tmp_path):
    base = tmp_path / "Simulations"
    path = base / TEMPLATE_FILENAME
    base.mkdir()
    path.write_text("- just a list\n", encoding="utf-8")
    record = base / "aa0"
    record.mkdir()
    with pytest.raises(ValueError, match="mapping"):
        rules.inherited_block(record, base)


# ---------------------------------------------------------------------------
# The package accessor: Experiment.bioschema
# ---------------------------------------------------------------------------


def test_an_experiment_resolves_its_block_and_keeps_its_metadata_as_stored(rules, tmp_path, monkeypatch):
    import fairmd.lipids.experiment as experiment

    base = tmp_path / "experiments"
    write_template(base, {"license": CC_BY})
    write_template(base / "OrderParameters", {"keywords": [BIOPHYSICS]})
    record = base / "OrderParameters" / "10.1039" / "B508190D" / "1"
    record.mkdir(parents=True)
    stored = {"name": "C-H bond order parameters of a POPC bilayer", "keywords": ["POPC"]}
    readme = {"MOLAR_FRACTIONS": {"POPC": 1.0}, "TEMPERATURE": 303, "bioschema_properties": stored}
    (record / "README.yaml").write_text(yaml.safe_dump(readme), encoding="utf-8")
    (record / "POPC_OrderParameters.json").write_text('{"M_G1_M M_G1H1_M": [[-0.1, 0.02]]}', encoding="utf-8")
    monkeypatch.setattr(experiment, "FMDL_EXP_PATH", str(base))

    exp = experiment.OPExperiment("10.1039/B508190D/1", str(record))
    assert exp.bioschema == {"name": stored["name"], "keywords": [BIOPHYSICS, "POPC"], "license": CC_BY}
    assert exp.metadata["bioschema_properties"] == stored


def test_an_experiment_without_a_block_has_no_resolved_one(rules, tmp_path):
    import fairmd.lipids.experiment as experiment

    record = tmp_path / "1"
    record.mkdir()
    (record / "README.yaml").write_text(yaml.safe_dump({"MOLAR_FRACTIONS": {"POPC": 1.0}}), encoding="utf-8")
    (record / "POPC_OrderParameters.json").write_text("{}", encoding="utf-8")
    assert experiment.OPExperiment("1", str(record)).bioschema is None
