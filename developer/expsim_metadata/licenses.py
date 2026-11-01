"""
Resolving a licence statement onto the SPDX list.

Registries state a licence as a URL, a name, or an SPDX id, and any of the three
may be missing or wrong. This module turns whatever was said into the pair the
Bioschemas block carries -- a resolvable URI and, where it is recognised, the
SPDX identifier -- and keeps the fetched SPDX list in one index so a run looks
it up rather than fetching it per record.
"""

import json
import re

from .constants import (
    ACCESS_RIGHTS_PREFIX,
    SPDX_LICENSES_URL,
)
from .helpers import fetch_json


def canonical_license_url(url: str) -> str:
    """Reduce a licence URL to a form comparable across sources.

    DataCite returns ``.../by/4.0``, ``.../by/4.0/`` and ``.../by/4.0/legalcode``
    for one licence while SPDX lists only the legalcode variant, so both sides
    pass through here before being compared.
    """
    text = str(url).strip().lower()
    text = re.sub(r"^https?://", "", text)
    text = re.sub(r"^www\.", "", text)
    text = re.sub(r"/(legalcode|deed)(\.[a-z-]{2,7})?$", "", text.rstrip("/"))
    return text.rstrip("/")


class SpdxIndex:
    """The SPDX licence list, indexed for lookup by URL or by name."""

    def __init__(self, payload):
        self.version = payload.get("licenseListVersion")
        self.by_url = {}
        self.by_name = {}
        self.by_id = {}
        for licence in payload.get("licenses") or []:
            self.by_id[licence["licenseId"]] = licence
            name = (licence.get("name") or "").strip().lower()
            if name:
                self.by_name.setdefault(name, licence)
            for see_also in licence.get("seeAlso") or []:
                self.by_url.setdefault(canonical_license_url(see_also), []).append(licence)

    def resolve(self, uri: str):
        """Return ``(licence, ambiguous)`` for a licence URI, or ``(None, False)``.

        Dozens of URLs in the SPDX list map to more than one identifier
        (GPL-2.0 vs GPL-2.0-only vs GPL-2.0-or-later, and this databank does hit
        that), so the choice is made deterministically: drop deprecated ids,
        drop legacy ``+`` ids, prefer the ``-only`` variant, then take the first
        alphabetically. A rerun cannot flip the answer.
        """
        candidates = self.by_url.get(canonical_license_url(uri))
        if not candidates:
            return None, False
        live = [c for c in candidates if not c.get("isDeprecatedLicenseId")] or candidates
        live = [c for c in live if not c["licenseId"].endswith("+")] or live
        ambiguous = len({c["licenseId"] for c in live}) > 1
        only = [c for c in live if c["licenseId"].endswith("-only")]
        return sorted(only or live, key=lambda c: c["licenseId"])[0], ambiguous

    def resolve_name(self, name):
        return self.by_name.get((name or "").strip().lower())


def load_spdx(cache_dir):
    """Fetch the SPDX licence list, caching it beside the DOI responses."""
    path = cache_dir / "spdx-licenses.json"
    if path.is_file():
        return SpdxIndex(json.loads(path.read_text(encoding="utf-8")))
    payload = fetch_json(SPDX_LICENSES_URL)
    if payload is None:
        return None
    cache_dir.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")
    return SpdxIndex(payload)


def license_block(licence, asserted_uri):
    block = {
        "spdx": licence["licenseId"],
        "name": licence.get("name") or None,
        "url": licence.get("reference") or None,
    }
    if asserted_uri:
        block["sameAs"] = asserted_uri
    return block


def compact_license(licence) -> dict:
    """Get licence block without its absent values.

    An unresolved licence keeps an explicit ``spdx: null`` beside the URI that
    failed to resolve, so the gap stays visible rather than looking like a
    licence nobody asserted; everything else absent is simply dropped.
    """
    compact = {key: value for key, value in licence.items() if value is not None}
    compact.setdefault("spdx", None)
    return compact


def resolve_license(rights_list, spdx):
    """Pick a licence from a DataCite ``rightsList``.

    Not ``rightsList[0]``: ``info:eu-repo/semantics/openAccess`` is an access
    status rather than a licence and is the most frequent entry in this
    databank, so taking the first element would stamp most records with a bogus
    licence. Take the first entry that actually resolves to an SPDX identifier.
    """
    fallback = None
    ambiguous = False
    for entry in rights_list or []:
        uri = (entry.get("rightsUri") or "").strip()
        name = (entry.get("rights") or "").strip()
        if uri.startswith(ACCESS_RIGHTS_PREFIX) or (not uri and not name):
            continue
        if spdx is not None:
            if uri:
                licence, ambiguous = spdx.resolve(uri)
                if licence:
                    return license_block(licence, uri), ambiguous
            licence = spdx.resolve_name(name)
            if licence:
                return license_block(licence, uri or None), False
        if fallback is None:
            # No SPDX match: the asserted URI is the only licence we have, so it
            # goes in `url` where schema.org/license expects it, with a null
            # `spdx` marking it as outside the vocabulary.
            fallback = {"spdx": None, "name": name or None, "url": uri or None}
    return fallback, ambiguous


def access_rights(rights_list: list) -> str | None:
    """Keep the COAR access status ``resolve_license`` skips, rather than lose it."""
    for entry in rights_list or []:
        uri = (entry.get("rightsUri") or "").strip()
        if uri.startswith(ACCESS_RIGHTS_PREFIX):
            return uri[len(ACCESS_RIGHTS_PREFIX):]
    return None
