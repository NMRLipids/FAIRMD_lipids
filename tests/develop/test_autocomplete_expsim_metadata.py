import importlib.util
import json
import sys
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft7Validator

# These tests exercise developer/autocomplete_expsim_metadata.py, which lives in
# the `developer/` folder and is not part of the distributed package. Mark the
# whole module as `develop` so it is isolated from the package test suite.
pytestmark = [pytest.mark.develop, pytest.mark.nodata]

ARTICLE_DOI = "10.1039/B508190D"
# The title the registry would answer with. Composed names must not contain it:
# every record of a hydration series would otherwise carry the same one.
ARTICLE_TITLE = "Structure and dynamics of DMPC bilayers"


def load_autocomplete_module():
    module_path = (
        Path(__file__).resolve().parents[2] / "developer" / "autocomplete_expsim_metadata.py"
    )
    spec = importlib.util.spec_from_file_location("autocomplete_expsim_metadata", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_databank(root, records):
    """A minimal databank tree: ``data_root_of`` needs both top-level folders."""
    (root / "Molecules" / "membrane" / "DMPC").mkdir(parents=True)
    (root / "Molecules" / "membrane" / "DMPC" / "metadata.yaml").write_text(
        yaml.safe_dump({"NMRlipids": {"name": "1,2-dimyristoyl-sn-glycero-3-phosphocholine"}}),
        encoding="utf-8",
    )
    (root / "Simulations").mkdir()

    paths = []
    for index, record in enumerate(records, start=1):
        directory = root / "experiments" / "OrderParameters" / "10.1039" / "B508190D" / str(index)
        directory.mkdir(parents=True)
        (directory / "DMPC_OrderParameters.json").write_text("{}", encoding="utf-8")
        path = directory / "README.yaml"
        path.write_text(yaml.safe_dump(record, sort_keys=False), encoding="utf-8")
        paths.append(path)
    return paths


def seed_cache(cache_dir):
    """A pre-baked CrossRef response, so the run touches no network."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "message": {
            "title": [ARTICLE_TITLE],
            "abstract": "<jats:p>An abstract that must not become a description.</jats:p>",
            "issued": {"date-parts": [[2005, 9]]},
            "publisher": "Royal Society of Chemistry",
            "container-title": ["Physical Chemistry Chemical Physics"],
            "author": [
                {"given": "Sergey V.", "family": "Dvinskikh"},
                {"given": "Oleg N.", "family": "Antzutkin"},
            ],
            "license": [],
        }
    }
    slug = ARTICLE_DOI.replace("/", "%2F")
    (cache_dir / f"{slug}.json").write_text(json.dumps([payload, "crossref"]), encoding="utf-8")


def base_record(**overrides):
    record = {
        "ARTICLE_DOI": ARTICLE_DOI,
        "TEMPERATURE": 314,
        "MEMBRANE_COMPOSITION": {"DMPC": 1},
        "TOTAL_HYDRATION": 95,
        "PH": "UNKNOWN",
        "NMR": {"METHOD": "2H:QE", "INSTRUMENT": "Bruker MSL-300"},
    }
    record.update(overrides)
    return record


@pytest.fixture
def generated(tmp_path):
    """Three sibling records that a fetched title could not tell apart."""
    mod = load_autocomplete_module()
    records = [
        base_record(),
        # Differs only in temperature.
        base_record(TEMPERATURE=324),
        # Differs only in the deprecated hydration field: no MEMBRANE_COMPOSITION,
        # no TOTAL_HYDRATION and no NMR block, which is how three real records in
        # BilayerData are written.
        {
            "ARTICLE_DOI": ARTICLE_DOI,
            "TEMPERATURE": 314,
            "MOLAR_FRACTIONS": {"DMPC": 1},
            "TOTAL_LIPID_CONCENTRATION": 11.1,
        },
    ]
    paths = build_databank(tmp_path, records)
    cache = tmp_path / "cache"
    seed_cache(cache)

    argv = ["autocomplete_expsim_metadata.py", "--cache", str(cache), *map(str, paths)]
    original_argv = sys.argv
    sys.argv = argv
    try:
        mod.main()
    finally:
        sys.argv = original_argv

    blocks = [
        yaml.safe_load(p.read_text(encoding="utf-8"))["bioschema_properties"] for p in paths
    ]
    return mod, paths, blocks


def test_names_are_composed_not_fetched(generated):
    _, _, blocks = generated
    for block in blocks:
        assert ARTICLE_TITLE not in block["name"]
        assert ARTICLE_TITLE not in block["description"]
        assert "abstract" not in block["description"].lower()
        assert block["name"].startswith("C-H bond order parameters of a DMPC bilayer")
        assert block["name"].endswith("[Dvinskikh 2005]")


def test_sibling_records_get_distinct_names(generated):
    _, _, blocks = generated
    names = [block["name"] for block in blocks]
    assert len(set(names)) == len(names)
    assert "at 314 K" in names[0]
    assert "at 324 K" in names[1]
    # The deprecated field is the only thing separating the third record.
    assert "5 waters per lipid" in names[2]


def test_legacy_composition_still_yields_keywords(generated):
    _, _, blocks = generated
    keywords = [k for k in blocks[2]["keywords"] if isinstance(k, str)]
    assert "DMPC" in keywords
    assert "1,2-dimyristoyl-sn-glycero-3-phosphocholine" in keywords


def test_fetched_title_moves_to_the_parent_article(generated):
    _, _, blocks = generated
    parent = blocks[0]["isPartOf"]
    assert parent["@type"] == "ScholarlyArticle"
    assert parent["name"] == ARTICLE_TITLE
    assert parent["isPartOf"]["name"] == "Physical Chemistry Chemical Physics"


def test_duplicate_names_are_reported(generated, tmp_path):
    mod, paths, _ = generated
    assert mod.duplicate_names(paths) == []

    # Force a clash and confirm the check catches it.
    clashed = yaml.safe_load(paths[1].read_text(encoding="utf-8"))
    clashed["bioschema_properties"]["name"] = yaml.safe_load(
        paths[0].read_text(encoding="utf-8")
    )["bioschema_properties"]["name"]
    paths[1].write_text(yaml.safe_dump(clashed, sort_keys=False), encoding="utf-8")

    found = mod.duplicate_names(paths)
    assert len(found) == 1
    assert found[0][0] == "ERROR"
    assert "duplicate name" in found[0][1]


def test_generated_block_is_schema_compliant(generated):
    _, paths, _ = generated
    # Read from the source tree rather than the installed package: this script is
    # developed against the schema in this checkout, and the test then runs
    # without the package being installed.
    schema_path = (
        Path(__file__).resolve().parents[2]
        / "src" / "fairmd" / "lipids" / "schema_validation" / "schema" / "experiment_schema.json"
    )
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    block_schema = schema["properties"]["bioschema_properties"]
    block_schema["definitions"] = schema["definitions"]

    for path in paths:
        block = yaml.safe_load(path.read_text(encoding="utf-8"))["bioschema_properties"]
        errors = sorted(Draft7Validator(block_schema).iter_errors(block), key=lambda e: e.path)
        assert not errors, f"{path}: {[e.message for e in errors]}"


def test_rerun_is_idempotent(generated):
    mod, paths, _ = generated
    before = [p.read_text(encoding="utf-8") for p in paths]

    root = paths[0].resolve().parents[5]
    names = mod.molecule_names(root)
    cache = root / "cache"
    spdx = mod.load_spdx(cache)
    for path in paths:
        assert mod.process(path, spdx, names, cache) is False

    assert [p.read_text(encoding="utf-8") for p in paths] == before
