"""
Script for membrane metadata autocomplete.

This script will try to fill further attributes in a
:ref:`membrane metadata file <addnewmol>` based on the information
queried with the inchikey from:
- UniChem
- ChEMBL
- ChEBI
- PubChem
- CAS Common Chemistry (requires the ``CAS_API_KEY`` environment variable)

.. note::
   This file is meant to be used by automated workflows.

   Several upstream services (notably EBI's UniChem and ChEBI) intermittently
   answer with transient ``5xx`` errors. Requests are therefore retried a few
   times with exponential backoff that honors any ``Retry-After`` header, so a
   blip does not abort metadata completion while staying polite to the servers.
   The retry budget can be overridden with the ``AUTOCOMPLETE_MAX_RETRIES``
   environment variable (set it to ``0`` to disable retries entirely).
"""

import json
import os
import random
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from html import unescape

import yaml

# HTTP status codes that signal a transient, server-side problem and are safe
# to retry. 429 (Too Many Requests) is included so we back off politely instead
# of hammering a rate-limited endpoint.
RETRYABLE_STATUS = frozenset({429, 500, 502, 503, 504})
MAX_RETRIES = max(0, int(os.environ.get("AUTOCOMPLETE_MAX_RETRIES", "4")))
BACKOFF_BASE = 1.0  # seconds for the first retry; doubles each attempt
MAX_BACKOFF = 30.0  # cap any single sleep so a flaky service can't stall us forever
DEFAULT_TIMEOUT = 15 
USER_AGENT = "FAIRMD-lipids-autocomplete (+https://github.com/NMRLipids/FAIRMD_lipids)"


def _retry_delay(error, attempt):
    """Seconds to wait before the next attempt.

    Prefers a server-provided ``Retry-After`` header (the polite signal), and
    otherwise falls back to exponential backoff with a little jitter so
    concurrent callers don't retry in lockstep.
    """
    headers = getattr(error, "headers", None)
    retry_after = headers.get("Retry-After") if headers is not None else None
    if retry_after:
        try:
            # Retry-After is usually a number of seconds; it may also be an HTTP
            # date, in which case we fall through to plain backoff.
            return min(float(retry_after), MAX_BACKOFF)
        except (TypeError, ValueError):
            pass
    backoff = BACKOFF_BASE * (2**attempt)
    return min(backoff, MAX_BACKOFF) + random.uniform(0, 0.5)


def fetch(req, timeout=DEFAULT_TIMEOUT):
    """Open ``req`` (a URL string or :class:`urllib.request.Request`) robustly.

    Returns the response body as ``bytes`` for an HTTP 200 response, or ``None``
    when the resource is unavailable. Transient failures (HTTP 429/5xx and
    connection-level errors such as timeouts) are retried with backed-off,
    ``Retry-After``-aware delays; definitive errors (e.g. 404) are not retried.
    """
    if isinstance(req, str):
        req = urllib.request.Request(req)
    req.add_header("User-Agent", USER_AGENT)

    for attempt in range(MAX_RETRIES + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read() if response.status == 200 else None
        except urllib.error.HTTPError as error:
            if error.code not in RETRYABLE_STATUS or attempt == MAX_RETRIES:
                return None
            delay = _retry_delay(error, attempt)
        except (urllib.error.URLError, TimeoutError):
            # Covers DNS failures, dropped connections and socket timeouts.
            if attempt == MAX_RETRIES:
                return None
            delay = _retry_delay(None, attempt)
        except Exception:
            return None
        time.sleep(delay)
    return None


def fetch_json(req, timeout=DEFAULT_TIMEOUT):
    """Like :func:`fetch`, but decode the body as JSON. Returns ``None`` on any
    failure (request error or malformed payload)."""
    body = fetch(req, timeout=timeout)
    if not body:
        return None
    try:
        return json.loads(body.decode("utf-8"))
    except (ValueError, UnicodeDecodeError):
        return None


def get_chembl(inchikey):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule?standard_inchi_key={inchikey}&format=json"
    return fetch_json(url) or {}


def get_pubchem(inchikey):
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
        f"{inchikey}/property/IUPACName,SMILES,InChI,InChIKey,MolecularFormula,MolecularWeight/JSON"
    )
    data = fetch_json(url)
    if data:
        try:
            return data["PropertyTable"]["Properties"][0]
        except (KeyError, IndexError, TypeError):
            pass
    return {}


def get_pubchem_synonyms(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
    data = fetch_json(url)
    if data:
        return data.get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])
    return []


def get_chebi(chebi_id):
    if not chebi_id:
        return {}

    url = f"https://www.ebi.ac.uk/chebi/backend/api/public/compound/{chebi_id}/?only_ontology_parents=false&only_ontology_children=false"

    return fetch_json(url) or {}


def get_metabolights(chebi_id):
    # MetaboLights reference compounds use the accession MTBLC<chebi numeric id>.
    # Return the identifier only if the compound exists in MetaboLights.
    if not chebi_id:
        return ""
    mtbl_id = f"MTBLC{chebi_id}"
    url = f"https://www.ebi.ac.uk/metabolights/ws/compounds/{mtbl_id}"
    return mtbl_id if fetch(url) is not None else ""


def get_cas(inchikey):
    # CAS Registry Numbers are not exposed by UniChem. They can be retrieved from
    # CAS Common Chemistry, which requires an API token supplied via the
    # CAS_API_KEY environment variable. Returns "" when the token is missing,
    # the service is unreachable, or no match is found.
    if not inchikey:
        return ""
    api_key = os.environ.get("CAS_API_KEY")
    if not api_key:
        return ""
    # CAS Common Chemistry requires field-qualified queries; a bare InChIKey
    # does not match, whereas "InChIKey=<value>" does.
    query = urllib.parse.quote(f"InChIKey={inchikey}")
    url = f"https://commonchemistry.cas.org/api/search?q={query}"
    req = urllib.request.Request(url, headers={"X-Api-Key": api_key})
    data = fetch_json(req)
    if data:
        results = data.get("results", [])
        if results:
            return results[0].get("rn", "") or ""
    return ""


def get_unichem(inchikey):
    url = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
    payload = json.dumps({"type": "inchikey", "compound": inchikey}).encode("utf-8")
    req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json"}, method="POST")
    data = fetch_json(req)
    if data:
        compounds = data.get("compounds", [])
        if compounds and "sources" in compounds[0]:
            return compounds[0]["sources"]
    return []


def extract_sameas(sources):
    mapping = {
        "pubchem": "pubchem.compound",
        "chebi": "ChEBI",
        "chembl": "ChEMBL",
        "lipidmaps": "lipidmaps",
        "metabolights": "metabolights",
        "swisslipids": "slm",
        "rcsb_pdb": "pdb.ligand",
        "pdbe": "pdb.ligand",
        "fdasrs": "unii",
        "cas": "cas",
    }
    result = {}
    for src in sources:
        prefix = mapping.get(src["shortName"])
        if prefix:
            value = src["compoundId"]
            if prefix == "ChEBI":
                if value:
                    value = value if str(value).startswith("CHEBI:") else f"CHEBI:{value}"
                else:
                    value = ""
            elif prefix == "pubchem.compound":
                try:
                    value = int(value)
                except ValueError:
                    pass
            result[prefix] = value
    return result


def get_chembl_id_from_unichem(sources):
    for src in sources:
        if src["shortName"] == "chembl":
            return src["compoundId"]
    return None


def clean_text(value):
    if not isinstance(value, str):
        return value
    return re.sub(r"<[^>]+>", "", unescape(value)).strip()


def safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def sanitize_sameas(sameas):
    patterns = {
        "ChEBI": r"^CHEBI:\d+$",
        "ChEMBL": r"^CHEMBL\d+$",
        "lipidmaps": r"^LM(FA|GL|GP|SP|ST|PR|SL|PK)[0-9]{4}([0-9a-zA-Z]{4,6})?$",
        "metabolights": r"^MTBL[CS]\d+$",
        "slm": r"^SLM:\d+$",
        "pdb.ligand": r"^[A-Za-z0-9]+$",
        "unii": r"^[A-Z0-9]+$",
        "cas": r"^\d{1,7}-\d{2}-\d$",
    }
    sanitized = {}
    for key, value in sameas.items():
        if key == "pubchem.compound":
            if isinstance(value, int):
                sanitized[key] = value
            else:
                try:
                    sanitized[key] = int(value)
                except (TypeError, ValueError):
                    print(
                        f"Warning: discarding sameAs '{key}' value {value!r}: not a valid integer.",
                        file=sys.stderr,
                    )
            continue
        pattern = patterns.get(key, r".+")
        if isinstance(value, str) and re.match(pattern, value):
            sanitized[key] = value
        else:
            print(
                f"Warning: discarding sameAs '{key}' value {value!r}: does not match expected pattern {pattern!r}.",
                file=sys.stderr,
            )
    return sanitized


def load_existing_metadata(path):
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)
    return {}


def update_metadata(existing, new_data):
    for key, value in new_data.items():
        if isinstance(value, dict):
            updated = update_metadata(existing.get(key, {}), value)
            if updated:
                existing[key] = updated
        elif isinstance(value, list):
            if value:
                existing[key] = existing.get(key, []) or value
        else:
            if value not in [None, "", {}]:
                existing[key] = existing.get(key) or value
    return existing


def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <metadata.yaml path>")
        sys.exit(1)

    metadata_path = sys.argv[1]

    # Extract NMRlipidsID from path (assumes structure: Molecules/membrane/<NMRlipidsID>/metadata.yaml)
    try:
        nmr_id = os.path.basename(os.path.dirname(metadata_path))
    except Exception:
        print("Error: Could not extract NMRlipidsID from path.")
        sys.exit(1)

    existing = load_existing_metadata(metadata_path)
    try:
        inchikey = existing["bioschema_properties"]["inChIKey"]
    except Exception:
        print("Error: Could not find bioschema_properties -> inChIKey in YAML file.")
        sys.exit(1)

    chembl = get_chembl(inchikey)
    pubchem = get_pubchem(inchikey)
    sources = get_unichem(inchikey)
    sameas = sanitize_sameas(extract_sameas(sources))

    cid = pubchem.get("CID", sameas.get("pubchem.compound"))
    synonyms = get_pubchem_synonyms(cid) if cid else []

    # First, check if there's a ChEBI ID from unichem
    chebi_id = sameas.get("ChEBI", "").replace("CHEBI:", "")
    chebi_data = get_chebi(chebi_id) if chebi_id else {}

    # MetaboLights is not exposed by UniChem; derive it from the ChEBI id.
    if chebi_id and "metabolights" not in sameas:
        metabolights_id = get_metabolights(chebi_id)
        if metabolights_id:
            sameas["metabolights"] = metabolights_id

    # CAS Registry Numbers are not exposed by UniChem; fetch them from CAS
    # Common Chemistry (requires the CAS_API_KEY environment variable).
    if "cas" not in sameas:
        cas_rn = get_cas(inchikey)
        if cas_rn and re.match(r"^\d{1,7}-\d{2}-\d$", cas_rn):
            sameas["cas"] = cas_rn

    # Collect alternate names with priority
    alternate_names = []

    # 1. Try ChEBI synonyms first
    if chebi_data and "names" in chebi_data:
        # Extract only the 'name' from SYNONYM type
        alternate_names = [
            syn["name"]
            for syn in chebi_data.get("names", {}).get("SYNONYM", [])
            if syn.get("type") == "SYNONYM" and syn.get("name")
        ]

    # 2. If no ChEBI synonyms, try ChEMBL synonyms
    if not alternate_names and chembl.get("molecule_synonyms"):
        alternate_names = [syn.get("molecule_synonym", "") for syn in chembl.get("molecule_synonyms", [])]

    # 3. If still no synonyms, try PubChem synonyms
    if not alternate_names and synonyms:
        alternate_names = synonyms
    alternate_names = [clean_text(name) for name in alternate_names if clean_text(name)]

    molecule_props = chembl.get("molecule_properties", {})
    molecule_structures = chembl.get("molecule_structures", {})
    chembl_id = get_chembl_id_from_unichem(sources)

    # Image selection logic
    if chembl_id:
        image_url = f"https://www.ebi.ac.uk/chembl/api/data/image/{chembl_id}?dimensions=200"
    elif cid:
        image_url = f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={cid}&t=l"
    else:
        image_url = ""

    nmr_name = (
        existing.get("NMRlipids", {}).get("name")
        or clean_text(chembl.get("pref_name", ""))
        or clean_text(molecule_props.get("iupac_name", ""))
        or clean_text(pubchem.get("IUPACName", ""))
        or nmr_id
    )

    bioschema = {
        "name": clean_text(molecule_props.get("iupac_name")) or clean_text(pubchem.get("IUPACName", "")),
        "iupacName": clean_text(molecule_props.get("iupac_name")) or clean_text(pubchem.get("IUPACName", "")),
        "molecularFormula": molecule_props.get("full_molformula") or pubchem.get("MolecularFormula", ""),
        "molecularWeight": safe_float(molecule_props.get("full_mwt") or pubchem.get("MolecularWeight")),
        "inChI": clean_text(molecule_structures.get("standard_inchi")) or clean_text(pubchem.get("InChI", "")),
        "inChIKey": clean_text(molecule_structures.get("standard_inchi_key"))
        or clean_text(pubchem.get("InChIKey", "")),
        "smiles": clean_text(molecule_structures.get("canonical_smiles")) or clean_text(pubchem.get("SMILES", "")),
        "image": image_url,
        "description": "",
    }

    if alternate_names:
        bioschema["alternateName"] = alternate_names

    new_data = {
        "NMRlipids": {"id": nmr_id, "name": nmr_name},
        "sameAs": sameas,
        "bioschema_properties": bioschema,
    }

    updated = update_metadata(existing, new_data)

    with open(metadata_path, "w", encoding="utf-8") as f:
        yaml.dump(updated, f, sort_keys=False, allow_unicode=True, default_flow_style=False)

    print(f"Updated metadata written to {metadata_path}")


if __name__ == "__main__":
    main()
