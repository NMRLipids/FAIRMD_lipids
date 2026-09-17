"""
Script for molecule metadata autocomplete.

This script will try to fill further attributes in a
:ref:`membrane metadata file <addnewmol>` based on the information
queried with the inchikey from:
- UniChem
- ChEMBL
- ChEBI
- PubChem
- CAS Common Chemistry (requires the ``CAS_API_KEY`` environment variable)
- LIPID MAPS (the structure depiction)

Structure depictions carry an ``imageAttribution`` block naming the source and its
licence, since all three image providers require a credit line. LIPID MAPS and
ChEMBL are credited under their Creative Commons licences; PubChem, which grants no
CC licence, is credited as its citation guidelines prescribe -- the sectioned record
URL, the CID and the primary PubChem paper. Pass ``--images-only`` to refresh just
the depiction and its credit from the identifiers already in the file, leaving every
other field and its formatting alone.

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

# Where structure depictions come from, in preference order, with everything the
# credit line needs: who published the picture, the record it belongs to, and the
# terms it is offered under. All three providers ask to be credited, so an image
# URL is never written without one.
#
# ``credit`` is the plain-text attribution. For the two Creative Commons sources it
# follows the CC recommended practices -- title, author, licence, source, every URL
# spelled out so the sentence stands on its own where markup cannot follow it:
# https://wiki.creativecommons.org/wiki/Recommended_practices_for_attribution
# PubChem is not under a CC licence and publishes its own citation guidelines
# instead, so its entry follows those rather than the CC phrasing.
#
# Optional keys, written only by the sources that define them:
#   ``depiction``  what the picture is called in the credit; PubChem asks for the
#                  "2D structure image" to be named as such.
#   ``identifier`` the source's own citable form of the id, e.g. ``CID 2244``.
#   ``usageInfo``  the terms-of-reuse page, when it is not the licence deed.
#   ``citation``   the paper the source asks to be cited alongside its data.
IMAGE_SOURCES = {
    "lipidmaps": {
        "name": "LIPID MAPS®",
        "url": "https://www.lipidmaps.org/",
        "image": "https://lipidmaps.org/api/molecules/{id}/svg",
        "record": "https://www.lipidmaps.org/databases/lmsd/{id}",
        "depiction": "structure depiction",
        "credit": (
            '"{title}" by LIPID MAPS® is licensed under CC BY 4.0 '
            "(https://creativecommons.org/licenses/by/4.0/). Source: {record}"
        ),
        "license": {
            "spdx": "CC-BY-4.0",
            "name": "Creative Commons Attribution 4.0 International",
            "url": "https://creativecommons.org/licenses/by/4.0/",
        },
    },
    "chembl": {
        "name": "ChEMBL",
        "url": "https://www.ebi.ac.uk/chembl/",
        "image": "https://www.ebi.ac.uk/chembl/api/data/image/{id}?dimensions=200",
        "record": "https://www.ebi.ac.uk/chembl/explore/compound/{id}",
        "depiction": "structure depiction",
        "credit": (
            '"{title}" by ChEMBL is licensed under CC BY-SA 3.0 '
            "(https://creativecommons.org/licenses/by-sa/3.0/). Source: {record}"
        ),
        "license": {
            "spdx": "CC-BY-SA-3.0",
            "name": "Creative Commons Attribution-ShareAlike 3.0 Unported",
            "url": "https://creativecommons.org/licenses/by-sa/3.0/",
        },
    },
    "pubchem": {
        "name": "PubChem",
        "url": "https://pubchem.ncbi.nlm.nih.gov/",
        "image": "https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={id}&t=l",
        # PubChem asks that a structure image be referenced as a section of the
        # compound record rather than as the record itself, hence the suffix:
        # https://pubchem.ncbi.nlm.nih.gov/docs/citation-guidelines
        # #section=Reusing-the-2D-or-3D-structure-image-of-a-compound-or-substance-record
        "record": "https://pubchem.ncbi.nlm.nih.gov/compound/{id}#section=2D-Structure",
        "identifier": "CID {id}",
        "depiction": "2D structure depiction",
        # Built from PubChem's citation guidelines rather than the CC phrasing: it
        # names the record the way PubChem asks (identifier plus sectioned URL),
        # states the permission it grants for structure images, and carries the
        # primary PubChem citation the guidelines ask for.
        "credit": (
            '"{title}" from PubChem {identifier} ({record}), National Center for '
            "Biotechnology Information, U.S. National Library of Medicine. "
            "PubChem 2D and 3D structure images may be reused without special "
            "permission; see https://pubchem.ncbi.nlm.nih.gov/docs/citation-guidelines. "
            "Cite: {citation[name]} {citation[identifier]}"
        ),
        # Not a licence grant but a web policy, hence no SPDX identifier.
        "license": {
            "name": "NLM Copyright and Privacy Policies",
            "url": "https://www.nlm.nih.gov/web_policies.html",
        },
        "usageInfo": (
            "https://pubchem.ncbi.nlm.nih.gov/docs/citation-guidelines"
            "#section=Reusing-the-2D-or-3D-structure-image-of-a-compound-or-substance-record"
        ),
        # The primary PubChem citation, in the AMA form the guidelines print.
        "citation": {
            "name": "Kim S, Chen J, Cheng T, et al. PubChem 2025 update. Nucleic Acids Res. 2025;53(D1):D1516-D1525.",
            "identifier": "doi:10.1093/nar/gkae1059",
            "url": "https://doi.org/10.1093/nar/gkae1059",
        },
    },
}


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


def lipidmaps_image(lm_id):
    """The LIPID MAPS SVG depiction for ``lm_id``, or ``""`` when it has none.

    An LMSD identifier does not guarantee a structure: LMGP12019AAA (TMCL) is a
    valid accession whose ``/svg`` endpoint answers 404, so the URL is only handed
    back once the service has actually served it.
    """
    if not lm_id:
        return ""
    url = IMAGE_SOURCES["lipidmaps"]["image"].format(id=lm_id)
    return url if fetch(url) is not None else ""


def image_attribution(source_key, identifier, nmr_id):
    """The credit that has to travel with a depiction from ``source_key``.

    Every source contributes ``creditText``, ``license`` and ``source``; the
    optional ``usageInfo``, ``source.identifier`` and ``citation`` appear only
    where the provider asks for them, as PubChem's citation guidelines do.
    """
    source = IMAGE_SOURCES[source_key]
    record = source["record"].format(id=identifier)
    title = f"{nmr_id} {source['depiction']}"
    citation = source.get("citation")
    # The source's own citable spelling of the id, e.g. ``CID 2244``. Empty for
    # sources that ask for nothing beyond the record URL.
    credited_id = source.get("identifier", "").format(id=identifier)

    attribution = {
        "creditText": source["credit"].format(
            title=title, record=record, identifier=credited_id, citation=citation or {}
        ),
        "license": dict(source["license"]),
    }
    if "usageInfo" in source:
        attribution["usageInfo"] = source["usageInfo"]
    attribution["source"] = {"name": source["name"], "url": source["url"], "sameAs": record}
    if credited_id:
        attribution["source"]["identifier"] = credited_id
    if citation:
        attribution["citation"] = dict(citation)
    return attribution


def select_image(sameas, chembl_id, cid, nmr_id):
    """Choose a structure depiction and its credit, as ``(url, attribution)``.

    LIPID MAPS comes first: it is the lipid-specific authority, it serves SVG
    rather than a raster fixed at one size, and CC BY 4.0 is the least demanding
    of the three sets of terms. ChEMBL and PubChem stand in where it has nothing.
    Returns ``("", None)`` when no source can depict the molecule.
    """
    lipidmaps_id = sameas.get("lipidmaps")
    url = lipidmaps_image(lipidmaps_id)
    if url:
        return url, image_attribution("lipidmaps", lipidmaps_id, nmr_id)

    for source_key, identifier in (("chembl", chembl_id), ("pubchem", cid)):
        if identifier:
            url = IMAGE_SOURCES[source_key]["image"].format(id=identifier)
            return url, image_attribution(source_key, identifier, nmr_id)

    return "", None


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


def apply_image(metadata, url, attribution):
    """Write the depiction and its credit, replacing whatever was there before.

    Unlike every other field these two are derived from the identifiers rather
    than accumulated, so they are assigned outright: a molecule that has gained a
    LIPID MAPS id must lose its PubChem picture, and a credit must never outlive
    the image it describes. A lookup that resolved nothing leaves the pair alone,
    since an unreachable service is not evidence that a depiction is gone.
    """
    if not url:
        return
    bioschema = metadata.setdefault("bioschema_properties", {})
    bioschema["image"] = url
    # Rebuilt rather than assigned into, so the credit sits next to the image it
    # describes: in a file written before attribution existed the new key would
    # otherwise land at the very end, far from the URL it belongs to.
    reordered = {}
    for key, value in bioschema.items():
        if key == "imageAttribution":
            continue
        reordered[key] = value
        if key == "image":
            reordered["imageAttribution"] = attribution
    metadata["bioschema_properties"] = reordered


def write_metadata(metadata_path, metadata):
    with open(metadata_path, "w", encoding="utf-8") as f:
        yaml.dump(metadata, f, sort_keys=False, allow_unicode=True, default_flow_style=False)
    print(f"Updated metadata written to {metadata_path}")


def parse_args(argv):
    images_only = False
    paths = []
    for arg in argv[1:]:
        if arg == "--images-only":
            images_only = True
        else:
            paths.append(arg)
    if len(paths) != 1:
        print("Usage: python autocomplete_mol_metadata.py [--images-only] <metadata.yaml path>")
        sys.exit(1)
    return images_only, paths[0]


def refresh_image(metadata_path, nmr_id, existing):
    """Re-derive only the depiction, from the identifiers already in the file.

    The ``sameAs`` block is all that picking an image needs, so this skips every
    registry lookup and touches no other field -- which keeps a re-run over the
    whole databank to one request per molecule and a diff that is only about
    images.
    """
    sameas = existing.get("sameAs") or {}
    url, attribution = select_image(
        sameas, sameas.get("ChEMBL"), sameas.get("pubchem.compound"), nmr_id
    )
    if not url:
        print(f"No image source for {nmr_id}; {metadata_path} left unchanged.", file=sys.stderr)
        return
    apply_image(existing, url, attribution)
    write_metadata(metadata_path, existing)


def main():
    images_only, metadata_path = parse_args(sys.argv)

    # Extract NMRlipidsID from path (assumes structure: Molecules/membrane/<NMRlipidsID>/metadata.yaml)
    try:
        nmr_id = os.path.basename(os.path.dirname(metadata_path))
    except Exception:
        print("Error: Could not extract NMRlipidsID from path.")
        sys.exit(1)

    existing = load_existing_metadata(metadata_path)

    if images_only:
        refresh_image(metadata_path, nmr_id, existing)
        return

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

    image_url, attribution = select_image(sameas, chembl_id, cid, nmr_id)

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
        "imageAttribution": attribution,
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
    apply_image(updated, image_url, attribution)

    write_metadata(metadata_path, updated)


if __name__ == "__main__":
    main()
