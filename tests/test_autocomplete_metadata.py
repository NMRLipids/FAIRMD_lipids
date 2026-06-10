import importlib.util
import json
import sys
from pathlib import Path

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
            {"shortName": "chebi", "compoundId": "1234"},
        ],
    )
    monkeypatch.setattr(mod, "get_pubchem_synonyms", lambda _: [])
    monkeypatch.setattr(
        mod,
        "get_chebi",
        lambda _: {"names": {"SYNONYM": [{"type": "SYNONYM", "name": "1-<em>OD&lt;/small&gt;-glucopyranoside"}]}},
    )

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
