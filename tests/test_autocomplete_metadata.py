import importlib.util
import json
import os
import re
import sys
from pathlib import Path

import pytest
import yaml
from jsonschema import Draft7Validator


def load_autocomplete_module():
    module_path = Path(__file__).resolve().parents[1] / "developer" / "autocomplete_metadata.py"
    spec = importlib.util.spec_from_file_location("autocomplete_metadata", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_autocomplete_output_is_schema_compliant(tmp_path, monkeypatch):
    mod = load_autocomplete_module()

    metadata_path = tmp_path / "Molecules" / "membrane" / "BOGUS" / "metadata.yaml"
    metadata_path.parent.mkdir(parents=True)
    metadata_path.write_text(
        yaml.safe_dump(
            {
                "NMRlipids": {"id": "BOGUS"},
                "bioschema_properties": {"inChIKey": "HEGSGKPQLMEBJL-RKQHYHRCSA-N"},
            }
        ),
        encoding="utf-8",
    )

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
    monkeypatch.setattr(
        mod,
        "get_unichem",
        lambda _: [
            {"shortName": "chembl", "compoundId": "CHEMBL446037"},
            {"shortName": "chebi", "compoundId": "CHEBI:1234"},
            {"shortName": "rcsb_pdb", "compoundId": "BOG"},
            {"shortName": "fdasrs", "compoundId": "V109WUT6RL"},
        ],
    )
    monkeypatch.setattr(mod, "get_pubchem_synonyms", lambda _: [])
    monkeypatch.setattr(
        mod,
        "get_chebi",
        lambda _: {"names": {"SYNONYM": [{"type": "SYNONYM", "name": "1-<em>OD&lt;/small&gt;-glucopyranoside"}]}},
    )
    monkeypatch.setattr(mod, "get_metabolights", lambda _: "MTBLC1234")
    monkeypatch.setattr(mod, "get_cas", lambda _: "29836-26-8")

    monkeypatch.setattr(sys, "argv", ["autocomplete_metadata.py", str(metadata_path)])
    mod.main()

    generated = yaml.safe_load(metadata_path.read_text(encoding="utf-8"))
    schema_path = (
        Path(__file__).resolve().parents[1]
        / "src"
        / "fairmd"
        / "lipids"
        / "schema_validation"
        / "schema"
        / "metadata_schema.json"
    )
    schema = json.loads(schema_path.read_text(encoding="utf-8"))

    errors = sorted(Draft7Validator(schema).iter_errors(generated), key=lambda e: e.path)
    assert not errors
    assert generated["NMRlipids"]["name"] == "(2R,3S)-name"
    assert generated["bioschema_properties"]["smiles"] == "CCCCCCCCOC@H]1[C@@HCO)O)O)O"
    assert generated["bioschema_properties"]["alternateName"] == ["1-OD-glucopyranoside"]
    assert generated["sameAs"]["ChEBI"] == "CHEBI:1234"
    assert generated["sameAs"]["pdb.ligand"] == "BOG"
    assert generated["sameAs"]["unii"] == "V109WUT6RL"
    assert generated["sameAs"]["metabolights"] == "MTBLC1234"
    assert generated["sameAs"]["cas"] == "29836-26-8"


@pytest.mark.network
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
