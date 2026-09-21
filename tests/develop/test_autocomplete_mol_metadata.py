import importlib.util
import json
import os
import re
import sys
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft7Validator

# These tests exercise developer/autocomplete_mol_metadata.py, which lives in the
# `developer/` folder and is not part of the distributed package. Mark the whole
# module as `develop` so it is isolated from the package test suite.
pytestmark = [pytest.mark.develop, pytest.mark.nodata]

LIPIDMAPS_SVG = "https://lipidmaps.org/api/molecules/LMGP01010005/svg"
CHEMBL_IMAGE = "https://www.ebi.ac.uk/chembl/api/data/image/CHEMBL446037?dimensions=200"
PUBCHEM_IMAGE = "https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=7906&t=l"


def load_autocomplete_module():
    module_path = Path(__file__).resolve().parents[2] / "developer" / "autocomplete_mol_metadata.py"
    spec = importlib.util.spec_from_file_location("autocomplete_mol_metadata", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_schema():
    # Read from the source tree rather than the installed package: this script is
    # developed against the schema in this checkout, and the test then runs
    # without the package being installed.
    schema_path = (
        Path(__file__).resolve().parents[2]
        / "src" / "fairmd" / "lipids" / "schema_validation" / "schema" / "metadata_schema.json"
    )
    return json.loads(schema_path.read_text(encoding="utf-8"))


def write_metadata(tmp_path, name, metadata):
    path = tmp_path / "Molecules" / "membrane" / name / "metadata.yaml"
    path.parent.mkdir(parents=True)
    path.write_text(yaml.safe_dump(metadata), encoding="utf-8")
    return path


def run(mod, path, *flags):
    argv = ["autocomplete_mol_metadata.py", *flags, str(path)]
    original_argv = sys.argv
    sys.argv = argv
    try:
        mod.main()
    finally:
        sys.argv = original_argv
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def stub_registries(mod, monkeypatch, unichem_sources, svg_available=False):
    """Answer every lookup from canned data, so no test reaches the network.

    ``svg_available`` decides what LIPID MAPS says about the depiction: `fetch` is
    the only path left open, and it is what `lipidmaps_image` probes with.
    """
    monkeypatch.setattr(mod, "get_chembl", lambda _: {"molecule_properties": {}, "molecule_structures": {}})
    monkeypatch.setattr(
        mod,
        "get_pubchem",
        lambda _: {
            "CID": 7906,
            "IUPACName": "(2R,3S)-name",
            "MolecularFormula": "C14H28O6",
            "MolecularWeight": 292.37,
            "InChI": "InChI=1S/...",
            "InChIKey": "HEGSGKPQLMEBJL-RKQHYHRCSA-N",
            "SMILES": "CCCCCCCCO<a>C@H]1[C@@H</a>CO)O)O)O",
        },
    )
    monkeypatch.setattr(mod, "get_unichem", lambda _: unichem_sources)
    monkeypatch.setattr(mod, "get_pubchem_synonyms", lambda _: [])
    monkeypatch.setattr(
        mod,
        "get_chebi",
        lambda _: {"names": {"SYNONYM": [{"type": "SYNONYM", "name": "1-<em>OD&lt;/small&gt;-glucopyranoside"}]}},
    )
    monkeypatch.setattr(mod, "get_metabolights", lambda _: "MTBLC1234")
    monkeypatch.setattr(mod, "get_cas", lambda _: "29836-26-8")
    monkeypatch.setattr(mod, "fetch", lambda *_, **__: b"<svg/>" if svg_available else None)


def test_autocomplete_output_is_schema_compliant(tmp_path, monkeypatch):
    mod = load_autocomplete_module()

    metadata_path = write_metadata(
        tmp_path,
        "BOGUS",
        {
            "NMRlipids": {"id": "BOGUS"},
            "bioschema_properties": {"inChIKey": "HEGSGKPQLMEBJL-RKQHYHRCSA-N"},
        },
    )

    stub_registries(
        mod,
        monkeypatch,
        [
            {"shortName": "chembl", "compoundId": "CHEMBL446037"},
            {"shortName": "chebi", "compoundId": "CHEBI:1234"},
            {"shortName": "rcsb_pdb", "compoundId": "BOG"},
            {"shortName": "fdasrs", "compoundId": "V109WUT6RL"},
        ],
    )

    generated = run(mod, metadata_path)

    errors = sorted(Draft7Validator(load_schema()).iter_errors(generated), key=lambda e: e.path)
    assert not errors
    assert generated["NMRlipids"]["name"] == "(2R,3S)-name"
    assert generated["bioschema_properties"]["smiles"] == "CCCCCCCCOC@H]1[C@@HCO)O)O)O"
    assert generated["bioschema_properties"]["alternateName"] == ["1-OD-glucopyranoside"]
    assert generated["sameAs"]["ChEBI"] == "CHEBI:1234"
    assert generated["sameAs"]["pdb.ligand"] == "BOG"
    assert generated["sameAs"]["unii"] == "V109WUT6RL"
    assert generated["sameAs"]["metabolights"] == "MTBLC1234"
    assert generated["sameAs"]["cas"] == "29836-26-8"

    # No LIPID MAPS id, so ChEMBL depicts it and its licence is the one credited.
    attribution = generated["bioschema_properties"]["imageAttribution"]
    assert generated["bioschema_properties"]["image"] == CHEMBL_IMAGE
    assert attribution["license"]["spdx"] == "CC-BY-SA-3.0"
    assert attribution["source"]["name"] == "ChEMBL"
    assert attribution["source"]["sameAs"] == "https://www.ebi.ac.uk/chembl/explore/compound/CHEMBL446037"
    assert "BOGUS structure depiction" in attribution["creditText"]
    assert "https://creativecommons.org/licenses/by-sa/3.0/" in attribution["creditText"]


def test_lipidmaps_is_preferred_and_replaces_an_existing_image(tmp_path, monkeypatch):
    """A LIPID MAPS depiction outranks the one already in the file."""
    mod = load_autocomplete_module()

    metadata_path = write_metadata(
        tmp_path,
        "POPC",
        {
            "NMRlipids": {"id": "POPC"},
            "bioschema_properties": {
                "inChIKey": "WTJKGGKOPKCXLL-VYOBOKEXSA-N",
                "image": PUBCHEM_IMAGE,
            },
        },
    )

    stub_registries(
        mod,
        monkeypatch,
        [
            {"shortName": "chembl", "compoundId": "CHEMBL446037"},
            {"shortName": "lipidmaps", "compoundId": "LMGP01010005"},
        ],
        svg_available=True,
    )

    generated = run(mod, metadata_path)
    bioschema = generated["bioschema_properties"]

    assert bioschema["image"] == LIPIDMAPS_SVG
    assert bioschema["imageAttribution"] == {
        "creditText": (
            '"POPC structure depiction" by LIPID MAPS® is licensed under CC BY 4.0 '
            "(https://creativecommons.org/licenses/by/4.0/). Source: "
            "https://www.lipidmaps.org/databases/lmsd/LMGP01010005"
        ),
        "license": {
            "spdx": "CC-BY-4.0",
            "name": "Creative Commons Attribution 4.0 International",
            "url": "https://creativecommons.org/licenses/by/4.0/",
        },
        "source": {
            "name": "LIPID MAPS®",
            "url": "https://www.lipidmaps.org/",
            "sameAs": "https://www.lipidmaps.org/databases/lmsd/LMGP01010005",
        },
    }
    # The CC-licensed sources ask for nothing beyond the deed, so the keys PubChem
    # needs are absent rather than empty.
    assert "usageInfo" not in bioschema["imageAttribution"]
    assert "citation" not in bioschema["imageAttribution"]
    assert not list(Draft7Validator(load_schema()).iter_errors(generated))


def test_lipidmaps_id_without_a_structure_falls_back(tmp_path, monkeypatch):
    """The real TMCL case: a valid LMSD id whose /svg endpoint answers 404."""
    mod = load_autocomplete_module()

    metadata_path = write_metadata(
        tmp_path,
        "TMCL",
        {
            "NMRlipids": {"id": "TMCL"},
            "bioschema_properties": {"inChIKey": "HEGSGKPQLMEBJL-RKQHYHRCSA-N"},
        },
    )

    stub_registries(
        mod,
        monkeypatch,
        [
            {"shortName": "lipidmaps", "compoundId": "LMGP12019AAA"},
            {"shortName": "pubchem", "compoundId": "7906"},
        ],
        svg_available=False,
    )

    generated = run(mod, metadata_path)
    bioschema = generated["bioschema_properties"]

    assert bioschema["image"] == PUBCHEM_IMAGE
    # PubChem grants no CC licence, so the credit follows its citation guidelines
    # instead: the 2D-Structure section of the record, the CID spelled the way
    # PubChem cites it, the reuse permission and the primary PubChem paper.
    assert bioschema["imageAttribution"] == {
        "creditText": (
            '"TMCL 2D structure depiction" from PubChem CID 7906 '
            "(https://pubchem.ncbi.nlm.nih.gov/compound/7906#section=2D-Structure), "
            "National Center for Biotechnology Information, U.S. National Library of Medicine. "
            "PubChem 2D and 3D structure images may be reused without special permission; "
            "see https://pubchem.ncbi.nlm.nih.gov/docs/citation-guidelines. "
            "Cite: Kim S, Chen J, Cheng T, et al. PubChem 2025 update. "
            "Nucleic Acids Res. 2025;53(D1):D1516-D1525. doi:10.1093/nar/gkae1059"
        ),
        # A web policy rather than a licence grant, hence no SPDX identifier.
        "license": {
            "name": "NLM Copyright and Privacy Policies",
            "url": "https://www.nlm.nih.gov/web_policies.html",
        },
        "usageInfo": (
            "https://pubchem.ncbi.nlm.nih.gov/docs/citation-guidelines"
            "#section=Reusing-the-2D-or-3D-structure-image-of-a-compound-or-substance-record"
        ),
        "source": {
            "name": "PubChem",
            "url": "https://pubchem.ncbi.nlm.nih.gov/",
            "sameAs": "https://pubchem.ncbi.nlm.nih.gov/compound/7906#section=2D-Structure",
            "identifier": "CID 7906",
        },
        "citation": {
            "name": "Kim S, Chen J, Cheng T, et al. PubChem 2025 update. Nucleic Acids Res. 2025;53(D1):D1516-D1525.",
            "identifier": "doi:10.1093/nar/gkae1059",
            "url": "https://doi.org/10.1093/nar/gkae1059",
        },
    }
    assert not list(Draft7Validator(load_schema()).iter_errors(generated))


def test_images_only_rewrites_the_pair_and_nothing_else(tmp_path, monkeypatch):
    """`--images-only` re-derives the depiction from `sameAs` alone."""
    mod = load_autocomplete_module()

    original = {
        "NMRlipids": {"id": "POPC", "name": "keep me"},
        "sameAs": {"lipidmaps": "LMGP01010005", "pubchem.compound": 5497103},
        "bioschema_properties": {
            "name": "keep me too",
            "molecularFormula": "C42H82NO8P",
            "image": PUBCHEM_IMAGE,
            "imageAttribution": {
                "creditText": "stale credit",
                "license": {"name": "NLM Copyright and Privacy Policies",
                            "url": "https://www.nlm.nih.gov/web_policies.html"},
                "source": {"name": "PubChem", "url": "https://pubchem.ncbi.nlm.nih.gov/"},
            },
        },
    }
    metadata_path = write_metadata(tmp_path, "POPC", original)

    # Any registry call would be a bug: --images-only must work from the file.
    for name in ("get_chembl", "get_pubchem", "get_unichem", "get_chebi", "get_cas"):
        monkeypatch.setattr(mod, name, lambda *_, called=name: pytest.fail(f"{called} must not be called"))
    monkeypatch.setattr(mod, "fetch", lambda *_, **__: b"<svg/>")

    generated = run(mod, metadata_path, "--images-only")
    bioschema = generated["bioschema_properties"]

    assert bioschema["image"] == LIPIDMAPS_SVG
    # The whole credit is replaced, not merged: nothing of PubChem's survives.
    assert bioschema["imageAttribution"]["license"]["spdx"] == "CC-BY-4.0"
    assert bioschema["imageAttribution"]["source"]["name"] == "LIPID MAPS®"

    generated["bioschema_properties"].pop("image")
    generated["bioschema_properties"].pop("imageAttribution")
    original["bioschema_properties"].pop("image")
    original["bioschema_properties"].pop("imageAttribution")
    assert generated == original


def test_unresolvable_image_leaves_the_existing_one_alone(tmp_path, monkeypatch):
    """An unreachable service is not evidence that a depiction is gone."""
    mod = load_autocomplete_module()

    attribution = {
        "creditText": "an existing credit",
        "license": {"spdx": "CC-BY-4.0", "name": "Creative Commons Attribution 4.0 International",
                    "url": "https://creativecommons.org/licenses/by/4.0/"},
        "source": {"name": "LIPID MAPS®", "url": "https://www.lipidmaps.org/"},
    }
    metadata_path = write_metadata(
        tmp_path,
        "POPC",
        {
            "NMRlipids": {"id": "POPC"},
            "sameAs": {"lipidmaps": "LMGP01010005"},
            "bioschema_properties": {
                "name": "POPC",
                "molecularFormula": "C42H82NO8P",
                "image": LIPIDMAPS_SVG,
                "imageAttribution": attribution,
            },
        },
    )

    # LIPID MAPS is down and no other identifier can stand in.
    monkeypatch.setattr(mod, "fetch", lambda *_, **__: None)

    generated = run(mod, metadata_path, "--images-only")

    assert generated["bioschema_properties"]["image"] == LIPIDMAPS_SVG
    assert generated["bioschema_properties"]["imageAttribution"] == attribution


def test_autocomplete_sameas_from_live_apis():
    """Live end-to-end check that the real APIs yield the expected cross references.

    Uses beta-octyl D-glucopyranoside (BOG). Skipped automatically when the
    external services are unreachable.
    """
    mod = load_autocomplete_module()

    inchikey = "HEGSGKPQLMEBJL-RKQHYHRCSA-N"

    sources = mod.get_unichem(inchikey)
    if not sources:
        pytest.skip("UniChem API unreachable; skipping live network test.")

    sameas = mod.sanitize_sameas(mod.extract_sameas(sources))
    chebi_id = sameas.get("ChEBI", "").replace("CHEBI:", "")
    if chebi_id and "metabolights" not in sameas:
        metabolights_id = mod.get_metabolights(chebi_id)
        if metabolights_id:
            sameas["metabolights"] = metabolights_id

    expected = {
        "ChEBI": "CHEBI:41128",
        "pubchem.compound": 62852,
        "metabolights": "MTBLC41128",
        "pdb.ligand": "BOG",
        "ChEMBL": "CHEMBL446037",
    }
    for key, value in expected.items():
        assert sameas.get(key) == value, f"{key}: expected {value!r}, got {sameas.get(key)!r}"

    # CAS Common Chemistry requires an API token; only verify when CAS_API_KEY is set.
    if os.environ.get("CAS_API_KEY"):
        cas_rn = mod.get_cas(inchikey)
        if cas_rn:
            assert re.match(r"^\d{1,7}-\d{2}-\d$", cas_rn), f"unexpected CAS format: {cas_rn!r}"
