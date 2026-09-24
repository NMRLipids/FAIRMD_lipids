import importlib
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
# A raw data deposition beside the article. Never looked up: the article is what
# the registries are asked about, so no cache entry is seeded for this one.
DATA_DOI = "10.18710/ETWNCU"
# The title the registry would answer with. Composed names must not contain it:
# every record of a hydration series would otherwise carry the same one.
ARTICLE_TITLE = "Structure and dynamics of DMPC bilayers"


DEVELOPER_DIR = Path(__file__).resolve().parents[2] / "developer"


def load_autocomplete_module():
    """The CLI script, loaded by path: developer/ is not an installed package."""
    module_path = DEVELOPER_DIR / "autocomplete_expsim_metadata.py"
    # The script imports the expsim_metadata package beside it the way it does when run
    # from that folder, so the folder has to be importable here.
    if str(DEVELOPER_DIR) not in sys.path:
        sys.path.insert(0, str(DEVELOPER_DIR))
    spec = importlib.util.spec_from_file_location("autocomplete_expsim_metadata", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_expsim_module(name):
    """One of the modules the script is split across, e.g. ``expsim_metadata.fields``."""
    if str(DEVELOPER_DIR) not in sys.path:
        sys.path.insert(0, str(DEVELOPER_DIR))
    return importlib.import_module(name)


SCHEMA_DIR = Path(__file__).resolve().parents[2] / "src" / "fairmd" / "lipids" / "schema_validation" / "schema"


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
    """Pre-baked registry responses, so the run touches no network."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    # The SPDX list is fetched once per run for the dataset licence; cache it too,
    # or every test in this module silently depends on spdx.org being reachable.
    (cache_dir / "spdx-licenses.json").write_text(
        json.dumps(
            {
                "licenseListVersion": "3.28.0",
                "licenses": [
                    {
                        "licenseId": "CC-BY-4.0",
                        "name": "Creative Commons Attribution 4.0 International",
                        "reference": "https://spdx.org/licenses/CC-BY-4.0.html",
                        "seeAlso": ["https://creativecommons.org/licenses/by/4.0/legalcode"],
                    },
                    # A second licence, so a fetched one can be told apart from
                    # the dataset licence this repository asserts.
                    {
                        "licenseId": "CC0-1.0",
                        "name": "Creative Commons Zero v1.0 Universal",
                        "reference": "https://spdx.org/licenses/CC0-1.0.html",
                        "seeAlso": ["https://creativecommons.org/publicdomain/zero/1.0/legalcode"],
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
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


def run_generator(tmp_path, records):
    """Enrich ``records`` as an experiment databank; return module, paths, blocks."""
    mod = load_autocomplete_module()
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

    blocks = [yaml.safe_load(p.read_text(encoding="utf-8"))["bioschema_properties"] for p in paths]
    return mod, paths, blocks


@pytest.fixture
def generated(tmp_path):
    """Three sibling records that a fetched title could not tell apart."""
    return run_generator(
        tmp_path,
        [
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
        ],
    )


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
    clashed["bioschema_properties"]["name"] = yaml.safe_load(paths[0].read_text(encoding="utf-8"))[
        "bioschema_properties"
    ]["name"]
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
    schema = block_schema("experiment_schema.json")

    for path in paths:
        block = yaml.safe_load(path.read_text(encoding="utf-8"))["bioschema_properties"]
        errors = sorted(Draft7Validator(schema).iter_errors(block), key=lambda e: e.path)
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


# ---------------------------------------------------------------------------
# Which DOI gets cited
#
# docs/src/schemas/experiment_metadata.md and experiment_schema.json both state
# the rule: with both DOIs given the deposition is what gets cited, and the
# article becomes the parent work. The cases below pin each row of that table.
# ---------------------------------------------------------------------------


@pytest.fixture
def both_dois(tmp_path):
    """One record carrying an ARTICLE_DOI and a DATA_DOI, the documented clash."""
    _, _, blocks = run_generator(tmp_path, [base_record(DATA_DOI=DATA_DOI)])
    return blocks[0]


def test_data_doi_is_what_gets_cited_when_both_are_given(both_dois):
    assert both_dois["citation"] == [DATA_DOI]


def test_data_doi_is_what_same_as_points_at_when_both_are_given(both_dois):
    assert both_dois["sameAs"] == f"https://doi.org/{DATA_DOI}"


def test_description_names_the_cited_doi_not_the_article(both_dois):
    assert f"https://doi.org/{DATA_DOI}" in both_dois["description"]
    assert ARTICLE_DOI not in both_dois["description"]


def test_article_stays_the_parent_work_when_both_are_given(both_dois):
    """The two rules are about different relations, so neither DOI is lost."""
    parent = both_dois["isPartOf"]
    assert parent["@type"] == "ScholarlyArticle"
    assert parent["identifier"] == ARTICLE_DOI
    # The article is still the record that was looked up: it carries the authors
    # and the journal that a raw data deposition does not.
    assert both_dois["_source"]["doi"] == ARTICLE_DOI
    assert parent["name"] == ARTICLE_TITLE
    assert [c["name"] for c in both_dois["creator"]] == [
        "Sergey V. Dvinskikh",
        "Oleg N. Antzutkin",
    ]


def test_article_only_record_cites_the_article(generated):
    """The other row of the table: with no DATA_DOI the article is the citation."""
    _, _, blocks = generated
    assert blocks[0]["citation"] == [ARTICLE_DOI]
    assert blocks[0]["sameAs"] == f"https://doi.org/{ARTICLE_DOI}"


def test_a_stale_article_citation_is_corrected_on_rerun(tmp_path):
    """Blocks written before the rule cited the article; a rerun must fix them."""
    stale = base_record(DATA_DOI=DATA_DOI)
    stale["bioschema_properties"] = {
        "name": "stale",
        "description": "stale",
        "citation": [ARTICLE_DOI, "10.1021/ja00000000"],
    }
    _, _, blocks = run_generator(tmp_path, [stale])

    assert ARTICLE_DOI not in blocks[0]["citation"]
    assert DATA_DOI in blocks[0]["citation"]
    # An unrelated citation is a related publication, which the rule says
    # nothing about, so it survives.
    assert "10.1021/ja00000000" in blocks[0]["citation"]


def test_check_reports_a_block_that_cites_the_wrong_doi(tmp_path):
    """--check sweeps a whole databank for records a rerun has to revisit."""
    mod, paths, _ = run_generator(tmp_path, [base_record(DATA_DOI=DATA_DOI)])
    spdx = mod.load_spdx(tmp_path / "cache")
    assert mod.check(paths[0], spdx) == []

    # Put the record back into the state every pre-rule block is in.
    readme = yaml.safe_load(paths[0].read_text(encoding="utf-8"))
    readme["bioschema_properties"]["citation"] = [ARTICLE_DOI]
    readme["bioschema_properties"]["sameAs"] = f"https://doi.org/{ARTICLE_DOI}"
    paths[0].write_text(yaml.safe_dump(readme, sort_keys=False), encoding="utf-8")

    messages = [message for level, message in mod.check(paths[0], spdx) if level == "ERROR"]
    assert any(DATA_DOI in m and "citation does not include" in m for m in messages)
    assert any(ARTICLE_DOI in m and "outranks" in m for m in messages)
    assert any("sameAs" in m for m in messages)


def test_data_ref_is_not_treated_as_a_doi(tmp_path):
    """DATA_REF is free text: it is never resolved, cited or linked to."""
    record = base_record(DATA_REF="nmrXiv sample S-1234, no DOI assigned")
    _, _, blocks = run_generator(tmp_path, [record])

    assert blocks[0]["citation"] == [ARTICLE_DOI]
    assert blocks[0]["sameAs"] == f"https://doi.org/{ARTICLE_DOI}"
    assert "nmrXiv" not in blocks[0]["description"]
    assert "S-1234" not in yaml.safe_dump(blocks[0])


def test_doi_roles_are_resolved_per_record():
    """record_dois() is the one place the precedence is decided."""
    mod = load_expsim_module("expsim_metadata.fields")

    both = mod.record_dois({"ARTICLE_DOI": ARTICLE_DOI, "DATA_DOI": DATA_DOI}, "experiments")
    assert (both.lookup, both.cited) == (ARTICLE_DOI, DATA_DOI)

    data_only = mod.record_dois({"DATA_DOI": DATA_DOI}, "experiments")
    assert (data_only.lookup, data_only.cited, data_only.article) == (DATA_DOI, DATA_DOI, None)

    # The deprecated experiment spelling of ARTICLE_DOI, still on 29 records: it
    # is read so those records keep their citation until they are renamed.
    legacy = mod.record_dois({"DOI": ARTICLE_DOI}, "experiments")
    assert (legacy.article, legacy.cited) == (ARTICLE_DOI, ARTICLE_DOI)

    # ... which also holds unpublished/<slug> values that are not DOIs at all.
    unpublished = mod.record_dois({"DOI": "unpublished/ferreira2023"}, "experiments")
    assert unpublished == (None, None, None, None)


def test_the_deprecated_doi_key_is_reported_not_rewritten():
    """DOI is deprecated for experiments, and named wherever it is read."""
    mod = load_expsim_module("expsim_metadata.fields")

    reported = mod.deprecated_keys({"DOI": ARTICLE_DOI}, "experiments")
    assert set(reported) == {"DOI"}
    assert "ARTICLE_DOI" in reported["DOI"] and "DATA_DOI" in reported["DOI"]

    # An unresolvable placeholder is still the deprecated key.
    assert set(mod.deprecated_keys({"DOI": "unpublished/ferreira2023"}, "experiments")) == {"DOI"}

    # Renamed records have nothing to report.
    assert mod.deprecated_keys({"ARTICLE_DOI": ARTICLE_DOI, "DATA_DOI": DATA_DOI}, "experiments") == {}

    # A simulation's DOI is the Zenodo deposition, which simulation_schema.json
    # declares: the same key is not deprecated there.
    assert mod.deprecated_keys({"DOI": SIMULATION_DOI}, "simulations") == {}


def test_a_record_still_using_doi_keeps_its_citation_and_is_flagged(tmp_path, capsys):
    """Enrichment does not rewrite the key, so --check keeps reporting it."""
    legacy = base_record()
    legacy["DOI"] = legacy.pop("ARTICLE_DOI")
    mod, paths, blocks = run_generator(tmp_path, [legacy])

    # Nothing is lost while the record waits to be renamed.
    assert blocks[0]["citation"] == [ARTICLE_DOI]
    assert "deprecated: DOI" in capsys.readouterr().out

    # The key is still in the record: renaming it is an edit to hand-written
    # content, made in the data repository rather than by this tool.
    assert "DOI" in yaml.safe_load(paths[0].read_text(encoding="utf-8"))

    spdx = mod.load_spdx(tmp_path / "cache")
    warnings = [message for level, message in mod.check(paths[0], spdx) if level == "WARNING"]
    assert any("DOI is deprecated" in m for m in warnings)


def test_dois_are_normalised_before_use():
    """A DOI given as a URL or with a doi: prefix still resolves and is cited bare."""
    mod = load_expsim_module("expsim_metadata.fields")
    prefixed = mod.record_dois(
        {"ARTICLE_DOI": f"doi:{ARTICLE_DOI}", "DATA_DOI": f"https://doi.org/{DATA_DOI}"},
        "experiments",
    )
    assert (prefixed.lookup, prefixed.cited) == (ARTICLE_DOI, DATA_DOI)


# ---------------------------------------------------------------------------
# The generated file, not just the generated block
#
# The block schema is mirrored in two files: experiment_schema.json describes an
# experiment README and simulation_schema.json a simulation one. The generator
# writes both kinds, so both contracts are checked against real output rather
# than against a hand-written fixture.
# ---------------------------------------------------------------------------

SIMULATION_DOI = "10.5281/zenodo.4040423"
DEPOSITION_TITLE = "Simulation trajectories of DMPC bilayers"


def build_simulation(root, record):
    """A minimal simulation databank tree, laid out as BilayerData does."""
    (root / "Molecules" / "membrane" / "DMPC").mkdir(parents=True)
    (root / "Molecules" / "membrane" / "DMPC" / "metadata.yaml").write_text(
        yaml.safe_dump({"NMRlipids": {"name": "1,2-dimyristoyl-sn-glycero-3-phosphocholine"}}),
        encoding="utf-8",
    )
    directory = root / "Simulations" / "aa0" / "bb1" / "cc2" / "dd3"
    directory.mkdir(parents=True)
    path = directory / "README.yaml"
    path.write_text(yaml.safe_dump(record, sort_keys=False), encoding="utf-8")
    return path


def seed_datacite_cache(cache_dir):
    cache_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "data": {
            "attributes": {
                "titles": [{"title": DEPOSITION_TITLE}],
                "publisher": "Zenodo",
                "publicationYear": 2020,
                "dates": [{"dateType": "Issued", "date": "2020-09-21"}],
                "version": "1.0",
                "creators": [
                    {
                        "name": "Ollila, O. H. Samuli",
                        "nameIdentifiers": [
                            {
                                "nameIdentifierScheme": "ORCID",
                                "nameIdentifier": "https://orcid.org/0000-0002-8135-3562",
                            }
                        ],
                    }
                ],
                "rightsList": [
                    {
                        "rights": "Creative Commons Attribution 4.0 International",
                        "rightsUri": "https://creativecommons.org/licenses/by/4.0/legalcode",
                    }
                ],
                "subjects": [{"subject": "lipid bilayer"}],
                "relatedIdentifiers": [],
            }
        }
    }
    slug = SIMULATION_DOI.replace("/", "%2F")
    (cache_dir / f"{slug}.json").write_text(json.dumps([payload, "datacite"]), encoding="utf-8")


SIMULATION_RECORD = {
    "DOI": SIMULATION_DOI,
    "TRJ": [["run.xtc"]],
    "TPR": [["run.tpr"]],
    "PREEQTIME": 100,
    "TIMELEFTOUT": 10,
    "SOFTWARE": "gromacs",
    "SYSTEM": "128DMPC_5000SOL_303K",
    "COMPOSITION": {
        "DMPC": {"NAME": "DMPC", "MAPPING": "mappingDMPCcharmm.yaml", "COUNT": [64, 64]},
    },
    "FF": "CHARMM36",
    "AUTHORS_CONTACT": "Ollila, O. H. Samuli",
    "ID": 566,
    "TRAJECTORY_SIZE": 3979765856,
    "TRJLENGTH": 200000.0,
    "NUMBER_OF_ATOMS": 40000,
    "DATEOFRUNNING": "2020-09-21",
    "TEMPERATURE": 303,
    "SOFTWARE_VERSION": "5.0.4",
}


def run_simulation(tmp_path, **overrides):
    """Enrich one simulation README; returns the module, its path and its content."""
    mod = load_autocomplete_module()
    path = build_simulation(tmp_path, {**SIMULATION_RECORD, **overrides})
    cache = tmp_path / "cache"
    seed_cache(cache)
    seed_datacite_cache(cache)

    argv = ["autocomplete_expsim_metadata.py", "--cache", str(cache), str(path)]
    original_argv = sys.argv
    sys.argv = argv
    try:
        mod.main()
    finally:
        sys.argv = original_argv

    return mod, path, yaml.safe_load(path.read_text(encoding="utf-8"))


@pytest.fixture
def generated_simulation(tmp_path):
    """One enriched simulation README, and the whole file it was written into."""
    return run_simulation(tmp_path)


def readme_schema():
    return json.loads((SCHEMA_DIR / "simulation_schema.json").read_text(encoding="utf-8"))


def test_generated_simulation_readme_validates(generated_simulation):
    """The whole generated file, against the schema a simulation README is held to."""
    _, path, readme = generated_simulation
    errors = sorted(Draft7Validator(readme_schema()).iter_errors(readme), key=lambda e: e.path)
    assert not errors, [f"{list(e.absolute_path)}: {e.message}" for e in errors]


def test_generated_simulation_block_is_populated(generated_simulation):
    """A file that validates but carries nothing would pass the test above."""
    _, _, readme = generated_simulation
    block = readme["bioschema_properties"]
    assert block["name"].startswith("Molecular dynamics trajectory of a DMPC bilayer")
    assert "NMRlipids simulation" not in block["name"]
    assert DEPOSITION_TITLE not in block["name"]
    assert block["alternateName"] == "128DMPC_5000SOL_303K"
    assert block["sameAs"] == f"https://doi.org/{SIMULATION_DOI}"
    assert block["license"]["spdx"] == "CC-BY-4.0"
    assert block["isPartOf"]["identifier"] == SIMULATION_DOI
    # The deposition title names the parent, not this record.
    assert block["isPartOf"]["name"] == DEPOSITION_TITLE


def test_duplicate_simulation_names_only_warn(generated_simulation):
    """Similar simulations may share a name; their ID only exists after merge."""
    mod, path, _ = generated_simulation
    twin = path.parent.parent / "twin" / "README.yaml"
    twin.parent.mkdir()
    twin.write_text(path.read_text(encoding="utf-8"), encoding="utf-8")

    found = mod.duplicate_names([path, twin])
    assert len(found) == 1
    assert found[0][0] == "WARNING"
    assert "duplicate name" in found[0][1]


def test_simulation_deposition_is_not_cited_as_a_publication(generated_simulation):
    """The deposition is this record, so it belongs in isPartOf, not in citation."""
    _, _, readme = generated_simulation
    assert SIMULATION_DOI not in (readme["bioschema_properties"].get("citation") or [])


def test_generated_simulation_rerun_is_idempotent(generated_simulation):
    mod, path, _ = generated_simulation
    before = path.read_text(encoding="utf-8")

    root = mod.data_root_of(path)
    cache = root / "cache"
    assert mod.process(path, mod.load_spdx(cache), mod.molecule_names(root), cache) is False
    assert path.read_text(encoding="utf-8") == before


def test_the_two_schemas_agree_on_the_bioschema_block():
    """Both schema files carry the block, and both say to keep them in sync.

    Without this, a property added to one file and forgotten in the other makes
    a block valid as an experiment and invalid as a simulation, or the reverse.
    """
    readme = readme_schema()
    experiment = json.loads((SCHEMA_DIR / "experiment_schema.json").read_text(encoding="utf-8"))
    from_readme = readme["properties"]["bioschema_properties"]
    from_experiment = experiment["properties"]["bioschema_properties"]

    # articleLicense is genuinely experiment-only: it holds the *article's*
    # licence, which only the experiment generator separates from the dataset's.
    experiment_only = {"articleLicense"}
    assert set(from_experiment["properties"]) - set(from_readme["properties"]) == experiment_only
    assert set(from_readme["properties"]) - set(from_experiment["properties"]) == set()
    for name in set(from_readme["properties"]):
        assert from_readme["properties"][name] == from_experiment["properties"][name], name
    for key in ("required", "additionalProperties", "patternProperties"):
        assert from_readme.get(key) == from_experiment.get(key), key

    # The definitions those properties $ref must match too, or the sync above is
    # only skin deep. `doi` is experiment-local and not referenced by the block.
    shared = set(readme["definitions"]) & set(experiment["definitions"])
    assert set(readme["definitions"]) - shared == set()
    assert set(experiment["definitions"]) - shared == {"doi"}
    for name in shared:
        assert readme["definitions"][name] == experiment["definitions"][name], name


# ---------------------------------------------------------------------------
# The documented examples
#
# Both schema pages show a filled-in bioschema_properties block. An example that
# no longer validates is worse than no example: it is what a contributor copies.
# ---------------------------------------------------------------------------

DOCS_DIR = Path(__file__).resolve().parents[2] / "docs" / "src" / "schemas"


def documented_block(page):
    """The ``bioschema_properties`` block from the fenced YAML example on a page."""
    text = (DOCS_DIR / page).read_text(encoding="utf-8")
    blocks = [
        part.split("```", 1)[0]
        for part in text.split("```yaml\n")[1:]
        if part.lstrip().startswith("bioschema_properties:")
    ]
    assert len(blocks) == 1, f"{page}: expected one bioschema_properties example, got {len(blocks)}"
    return yaml.safe_load(blocks[0])["bioschema_properties"]


def block_schema(filename):
    """The bioschema_properties subschema, made standalone so it can be validated."""
    schema = json.loads((SCHEMA_DIR / filename).read_text(encoding="utf-8"))
    block = schema["properties"]["bioschema_properties"]
    block["definitions"] = schema["definitions"]
    return block


@pytest.mark.parametrize(
    ("page", "schema_file"),
    [
        ("experiment_metadata.md", "experiment_schema.json"),
        ("simulation_metadata.md", "simulation_schema.json"),
    ],
)
def test_documented_example_validates(page, schema_file):
    block = documented_block(page)
    errors = sorted(Draft7Validator(block_schema(schema_file)).iter_errors(block), key=lambda e: e.path)
    assert not errors, f"{page}: {[f'{list(e.absolute_path)}: {e.message}' for e in errors]}"


@pytest.mark.parametrize(
    ("page", "ending"),
    [
        # Experiments end in their source tag, simulations in their method detail.
        ("experiment_metadata.md", "]"),
        ("simulation_metadata.md", ")"),
    ],
)
def test_documented_example_shows_a_composed_name(page, ending):
    """The examples must not go back to showing a fetched registry title."""
    block = documented_block(page)
    assert block["name"].rstrip().endswith(ending), block["name"]
    # The composed description ends by naming where the record came from.
    assert block["description"].rstrip().endswith("."), block["description"]


# ---------------------------------------------------------------------------
# What a run must not destroy
#
# A registry that answered yesterday and times out today must cost a record
# nothing: the block it already carries is the only copy of the creators, dates
# and licence that were fetched, and a rewritten block cannot get them back.
# ---------------------------------------------------------------------------


DEPOSITION_ONLY_DOI = "10.57992/nmrxiv.p157.s1600"


def seed_deposition_cache(cache_dir, doi=DEPOSITION_ONLY_DOI):
    """An nmrXiv-style deposition: DataCite answers, and there is no article."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "data": {
            "attributes": {
                "titles": [{"title": "Raw NMR data for a cardiolipin bilayer"}],
                "publisher": "nmrXiv",
                "publicationYear": 2023,
                "dates": [{"dateType": "Issued", "date": "2023-05-04"}],
                "creators": [{"name": "Hörstel, Marie"}],
                # Deliberately not the licence this repository distributes under,
                # so the two cannot be confused in the assertions below.
                "rightsList": [
                    {
                        "rights": "Creative Commons Zero v1.0 Universal",
                        "rightsUri": "https://creativecommons.org/publicdomain/zero/1.0/legalcode",
                    }
                ],
                "subjects": [],
                "relatedIdentifiers": [],
            }
        }
    }
    slug = doi.replace("/", "%2F")
    (cache_dir / f"{slug}.json").write_text(json.dumps([payload, "datacite"]), encoding="utf-8")


def test_a_failed_lookup_leaves_the_last_good_block_alone(generated, tmp_path, monkeypatch):
    """A transient outage must not strip an already enriched record."""
    mod, paths, _ = generated
    registries = load_expsim_module("expsim_metadata.registries")
    before = [p.read_text(encoding="utf-8") for p in paths]

    root = paths[0].resolve().parents[5]
    names = mod.molecule_names(root)
    spdx = mod.load_spdx(root / "cache")

    # An empty cache and a registry that answers nothing: the lookup fails the
    # way an outage makes it fail.
    monkeypatch.setattr(registries, "fetch_json", lambda *args, **kwargs: None)
    empty_cache = tmp_path / "empty-cache"
    for path in paths:
        assert mod.process(path, spdx, names, empty_cache) is False

    assert [p.read_text(encoding="utf-8") for p in paths] == before


def test_a_failed_lookup_is_not_cached(tmp_path, monkeypatch):
    """Caching a failure would suppress every later attempt at that DOI."""
    registries = load_expsim_module("expsim_metadata.registries")
    cache = tmp_path / "cache"

    monkeypatch.setattr(registries, "fetch_json", lambda *args, **kwargs: None)
    assert registries.resolve_doi(ARTICLE_DOI, "experiments", cache)[0] is None
    assert list(cache.glob("*.json")) == []

    payload = {"message": {"title": [ARTICLE_TITLE]}}
    monkeypatch.setattr(registries, "fetch_json", lambda *args, **kwargs: payload)
    assert registries.resolve_doi(ARTICLE_DOI, "experiments", cache)[0] == payload
    assert len(list(cache.glob("*.json"))) == 1


# ---------------------------------------------------------------------------
# Which licence belongs to what
# ---------------------------------------------------------------------------


def test_a_deposition_licence_is_not_called_an_article_licence(tmp_path):
    """With only a DATA_DOI there is no article, so nothing may claim one.

    The four nmrXiv records in BilayerData are exactly this case: DataCite
    answers for a deposition, and its licence is a property of that deposition
    rather than of a paper that was never written.
    """
    mod = load_autocomplete_module()
    paths = build_databank(tmp_path, [base_record(DATA_DOI=DEPOSITION_ONLY_DOI, ARTICLE_DOI=None)])
    record = yaml.safe_load(paths[0].read_text(encoding="utf-8"))
    del record["ARTICLE_DOI"]
    paths[0].write_text(yaml.safe_dump(record, sort_keys=False), encoding="utf-8")

    cache = tmp_path / "cache"
    seed_cache(cache)
    seed_deposition_cache(cache)
    argv = ["autocomplete_expsim_metadata.py", "--cache", str(cache), str(paths[0])]
    original_argv = sys.argv
    sys.argv = argv
    try:
        mod.main()
    finally:
        sys.argv = original_argv

    block = yaml.safe_load(paths[0].read_text(encoding="utf-8"))["bioschema_properties"]
    assert "articleLicense" not in block
    # What this repository distributes stays the dataset licence ...
    assert block["license"]["spdx"] == "CC-BY-4.0"
    # ... and the deposition keeps its own, on the deposition.
    parent = block["isPartOf"]
    assert parent["identifier"] == DEPOSITION_ONLY_DOI
    assert parent["license"]["spdx"] == "CC0-1.0"

    errors = sorted(Draft7Validator(block_schema("experiment_schema.json")).iter_errors(block), key=lambda e: e.path)
    assert not errors, [e.message for e in errors]


# ---------------------------------------------------------------------------
# Retiring PUBLICATION
# ---------------------------------------------------------------------------


def test_publication_is_retired_for_simulations(tmp_path):
    """A simulation retires the field as well, and still validates without it.

    ``simulation_schema.json`` declares ``PUBLICATION`` but does not require it,
    so a record that has given it up is still a valid simulation. Keeping it
    would leave a superseded field alive beside the ``citation`` list that
    replaced it, free to drift out of step with it.
    """
    _, path, readme = run_simulation(tmp_path, PUBLICATION=f"https://doi.org/{ARTICLE_DOI}")
    assert "PUBLICATION" not in path.read_text(encoding="utf-8")
    assert ARTICLE_DOI in readme["bioschema_properties"]["citation"]
    errors = sorted(Draft7Validator(readme_schema()).iter_errors(readme), key=lambda e: e.path)
    assert not errors, [f"{list(e.absolute_path)}: {e.message}" for e in errors]


def test_publication_is_retired_for_experiments(tmp_path):
    """experiment_schema.json does not declare PUBLICATION and forbids extras."""
    _, paths, blocks = run_generator(
        tmp_path, [base_record(PUBLICATION=f"Dvinskikh et al., https://doi.org/{ARTICLE_DOI}")]
    )
    assert blocks[0]["citation"] == [ARTICLE_DOI]
    assert "PUBLICATION" not in paths[0].read_text(encoding="utf-8")


def test_free_text_in_publication_moves_to_citation_before_the_field_goes(tmp_path):
    """A reference with no DOI is kept verbatim, so nothing is lost with the field."""
    _, paths, blocks = run_generator(tmp_path, [base_record(PUBLICATION="Dvinskikh et al., PCCP 7 (2005) 3255")])
    assert "Dvinskikh et al., PCCP 7 (2005) 3255" in blocks[0]["citation"]
    assert "PUBLICATION" not in paths[0].read_text(encoding="utf-8")


def test_a_publication_the_citation_rule_removed_is_kept(tmp_path):
    """The field only goes once ``citation`` demonstrably carries its content.

    With both DOIs given the article is deliberately *not* cited -- it is the
    parent work instead -- so a PUBLICATION naming it is not represented there
    and stays where it is rather than being dropped on the strength of a rule
    that removed it.
    """
    _, paths, blocks = run_generator(
        tmp_path, [base_record(DATA_DOI=DATA_DOI, PUBLICATION=f"https://doi.org/{ARTICLE_DOI}")]
    )
    assert blocks[0]["citation"] == [DATA_DOI]
    assert "PUBLICATION" in paths[0].read_text(encoding="utf-8")


def test_a_block_scalar_publication_is_removed_whole(tmp_path):
    """A string may be written over several lines, and all of them have to go.

    Removing the ``PUBLICATION:`` line alone strands its indented text at top
    level, which either fails to parse or is swallowed by the next key.
    """
    records = load_expsim_module("expsim_metadata.records")
    readme = (
        "TEMPERATURE: 314\n"
        "PUBLICATION: |\n"
        f"  Dvinskikh et al., PCCP 7 (2005) 3255.\n"
        f"  https://doi.org/{ARTICLE_DOI}\n"
        "\n"
        "TOTAL_HYDRATION: 95\n"
    )
    out, removed = records.drop_publication(readme, {"citation": [ARTICLE_DOI]})

    assert removed
    parsed = yaml.safe_load(out)
    assert "PUBLICATION" not in parsed
    # The stranded lines did not attach themselves to the surviving keys.
    assert parsed == {"TEMPERATURE": 314, "TOTAL_HYDRATION": 95}


def test_a_multi_line_publication_the_rule_kept_is_left_intact(tmp_path):
    """Nothing is touched when the citation does not cover the field."""
    records = load_expsim_module("expsim_metadata.records")
    readme = "PUBLICATION: |\n  Dvinskikh et al., PCCP 7 (2005) 3255\nTOTAL_HYDRATION: 95\n"
    out, removed = records.drop_publication(readme, {"citation": [DATA_DOI]})

    assert not removed
    assert out == readme


# ---------------------------------------------------------------------------
# Where a record sits
# ---------------------------------------------------------------------------


def test_the_toy_data_simulation_folders_count_as_simulations(tmp_path):
    """ToyData splits trajectories across Simulations.1, .2 and .AddData."""
    fields = load_expsim_module("expsim_metadata.fields")
    for folder in ("Simulations", "Simulations.1", "Simulations.2", "Simulations.AddData"):
        path = tmp_path / folder / "aa0" / "README.yaml"
        assert fields.record_kind(path) == "simulations", folder
    assert fields.record_kind(tmp_path / "experiments" / "OrderParameters" / "1" / "README.yaml") == "experiments"


def test_the_databank_root_is_found_in_the_toy_layout(tmp_path):
    fields = load_expsim_module("expsim_metadata.fields")
    root = tmp_path / "ToyData"
    (root / "Molecules" / "membrane").mkdir(parents=True)
    (root / "Simulations.1" / "aa0").mkdir(parents=True)
    assert fields.data_root_of(root / "Simulations.1" / "aa0" / "README.yaml") == root.resolve()


# ---------------------------------------------------------------------------
# The pH scale
# ---------------------------------------------------------------------------


def test_a_ph_range_is_bounded_to_the_ph_scale():
    """Both endpoints obey the bound the numeric branch states, decimals included."""
    schema = json.loads((SCHEMA_DIR / "experiment_schema.json").read_text(encoding="utf-8"))
    validator = Draft7Validator(schema["properties"]["PH"])

    for value in (7, 7.4, 0, 14, "UNKNOWN", "8-10", "1.3-13.2", "0-14", "14.0-14.0", "3.0-9.25"):
        assert validator.is_valid(value), value
    for value in ("14.5-14.5", "15-16", "-1-5", "10-", "7,4", "8 - 10"):
        assert not validator.is_valid(value), value
