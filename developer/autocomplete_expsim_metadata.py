"""
Script for simulation and experiment Bioschemas autocomplete.

This script fills a ``bioschema_properties`` block in a simulation or experiment
``README.yaml`` so the record can be published as a `Bioschemas Dataset
<https://bioschemas.org/profiles/Dataset/1.0-RELEASE>`_. It is the dataset-level
counterpart of ``autocomplete_mol_metadata.py``, which does the same job for
molecules.

``name`` and ``description`` are composed from the record's own fields --
composition, temperature, hydration, ions, method -- and never taken from the
registry. A registry title names the paper or the Zenodo deposition, not the one
measurement or trajectory the record holds, so it describes the wrong thing and
does not tell sibling records apart: before this, 104 experiments shared 37
titles and 896 simulations shared 479. The fetched title is kept where it is
true, as ``isPartOf.name`` on the parent work, and fetched abstracts are dropped.

Every composed title ends in a bracketed tag that keeps near-identical records
apart -- the databank ``ID`` for a simulation, the first author and year for an
experiment. ``--check`` reports any two records that still share a title, and is
meant to be pointed at the whole databank, since only a full run can see a clash.

An experiment can carry two DOIs that mean different things, and ``record_dois``
decides once which is used for what: the ``ARTICLE_DOI`` is looked up, because
CrossRef holds the authors and the journal, and it becomes ``isPartOf``; the
``DATA_DOI`` is what gets cited and what ``sameAs`` points at, because the rule
in ``docs/src/schemas/experiment_metadata.md`` is to cite the data. With only one
of them given, that one fills every role.

The remaining values are resolved from the record's DOI:

- DataCite -- Zenodo depositions, i.e. every simulation
- CrossRef -- journal articles, i.e. most experiments
- SPDX     -- the licence list, used as the licence controlled vocabulary

giving creators, dates, licence, publisher, citations and the parent work's
title. These are enriched from data already in the repository (composition,
force field, NMR/X-ray method, the analysis outputs present beside the README)
with terms from EDAM, CHMO and UO.

Properties that depend on where the databank is deployed -- ``identifier``,
``url``, ``@id``, ``@type``, ``@context``, ``dct:conformsTo`` and
``includedInDataCatalog`` -- are deliberately *not* written here; the web
frontend supplies them.

The script also normalises the surrounding record: dates are rewritten as quoted
ISO strings, and ``PUBLICATION`` is retired in favour of ``citation`` once its
content is safely represented there.

.. note::
   This file is meant to be used by automated workflows.

   Unlike ``autocomplete_mol_metadata.py`` the file is **not** re-serialised. The
   generated block is appended, and only that block is rewritten in place, so
   hand-written comments elsewhere in the README survive -- most experiment
   files carry them. Runs are idempotent: a record whose content has not changed
   keeps its original ``_source.retrieved`` date.

   Upstream services answer with transient ``5xx`` errors from time to time, so
   requests are retried with exponential backoff that honours any
   ``Retry-After`` header. The retry budget can be overridden with the
   ``AUTOCOMPLETE_MAX_RETRIES`` environment variable (``0`` disables retries).
"""

import argparse
import json
import os
import random
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import namedtuple
from datetime import date
from pathlib import Path

import yaml

# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------

RETRYABLE_STATUS = frozenset({429, 500, 502, 503, 504})
MAX_RETRIES = max(0, int(os.environ.get("AUTOCOMPLETE_MAX_RETRIES", "4")))
BACKOFF_BASE = 1.0
MAX_BACKOFF = 30.0
DEFAULT_TIMEOUT = 30
# The project contact address puts CrossRef requests in its polite pool.
USER_AGENT = (
    "FAIRMD-lipids-bioschema (+https://github.com/NMRLipids/FAIRMD_lipids; "
    "mailto:databank@nmrlipids.fi)"
)

DATACITE_URL = "https://api.datacite.org/dois/{doi}"
CROSSREF_URL = "https://api.crossref.org/works/{doi}"
SPDX_LICENSES_URL = "https://spdx.org/licenses/licenses.json"

# ---------------------------------------------------------------------------
# Vocabularies
# ---------------------------------------------------------------------------

EDAM_SET = "http://edamontology.org"
OBO_SET = "http://purl.obolibrary.org/obo"
QUDT_UNIT = "http://qudt.org/vocab/unit"

EDAM_SIMULATION = [
    ("topic_3892", "Biomolecular simulation"),
    ("topic_0176", "Molecular dynamics"),
    ("topic_0820", "Membrane and lipoproteins"),
    ("topic_3306", "Biophysics"),
    ("data_3870", "Trajectory data"),
]
EDAM_NMR = [
    ("topic_0593", "NMR"),
    ("topic_0820", "Membrane and lipoproteins"),
    ("topic_3306", "Biophysics"),
]
# EDAM has no small-angle-scattering term; X-ray diffraction is the closest.
EDAM_XRAY = [
    ("topic_2828", "X-ray diffraction"),
    ("topic_0820", "Membrane and lipoproteins"),
    ("topic_3306", "Biophysics"),
]

# CHMO is exact for the experimental methods and EDAM sits alongside it; neither
# covers this databank alone, since EDAM has no scattering term and CHMO no
# molecular-dynamics term. EDAM-Bioimaging was evaluated as a fallback and is
# not citable: it declares itself 1.0alpha_pre-pre-release and only 36 of its
# 505 classes carry stable numeric ids, those being EDAM terms it re-uses.
CHMO_SAXS = ("CHMO_0000204", "small-angle X-ray scattering")
CHMO_NMR = ("CHMO_0000591", "nuclear magnetic resonance spectroscopy")
CHMO_SSNMR = ("CHMO_0000614", "solid-state nuclear magnetic resonance spectroscopy")
CHMO_PDLF = ("CHMO_0001067", "proton detected separated local field spectroscopy")
CHMO_PDLF_R = ("CHMO_0001068", "R-type recoupling proton detected separated local field spectroscopy")

EDAM_MD_OP = ("operation_2476", "Molecular dynamics")
EDAM_NMR_TOPIC = ("topic_0593", "NMR")
EDAM_XRAY_TOPIC = ("topic_2828", "X-ray diffraction")

# Units. UO is preferred -- it is OBO, resolves through the same OLS endpoint as
# EDAM and CHMO, and has a square-angstrom term QUDT lacks -- but it has no
# reciprocal units, so the form factor's q axis falls back to QUDT.
#
# The density profile is a number density in nm^-3 and stays symbol-only on
# purpose: UO has nothing suitable (UO_0000177 is a solution-concentration
# grouping, UO_0000182 is mass density) and QUDT has no PER-NanoM3. QUDT's
# PER-M3 has the right dimension but is wrong by 10^27 at this scale, and a
# precise-looking wrong code is worse than an honest plain-text symbol.
UNITS = {
    "area per lipid": (f"{OBO_SET}/UO_0000324", "Å²"),
    "membrane thickness": (f"{OBO_SET}/UO_0000018", "nm"),
    "C-H bond order parameter": (f"{OBO_SET}/UO_0000186", "dimensionless"),
    "equilibration time": (f"{OBO_SET}/UO_0000150", "ns"),
    "X-ray scattering form factor": (f"{QUDT_UNIT}/PER-ANGSTROM", "Å⁻¹"),
    "total density profile": (None, "nm⁻³"),
}

# Licence of the data this repository distributes, per its own LICENSE file. An
# experiment record digitises values out of a paper: the paper carries the
# publisher's copyright, the digitised values are distributed under CC-BY-4.0.
# Those are different facts, so `license` states the dataset licence and
# `articleLicense` keeps whatever the registry said about the publication.
DATASET_LICENSE_SPDX = "CC-BY-4.0"
DATASET_LICENSE_URI = "https://creativecommons.org/licenses/by/4.0/"

# COAR access statuses. DataCite returns these in rightsList alongside real
# licences -- and they are the most frequent entry there -- but they describe
# access, not terms of reuse.
ACCESS_RIGHTS_PREFIX = "info:eu-repo/semantics/"

# DataCite relation types pointing at the paper describing a deposition.
# IsVersionOf is excluded: it points at Zenodo's own concept DOI, which is
# another copy of this dataset rather than a citation.
CITATION_RELATIONS = frozenset({"IsSupplementTo", "IsDocumentedBy", "IsCitedBy"})

# Which analysis outputs correspond to which measured quantity.
VARIABLES = {
    "apl.json": "area per lipid",
    "thickness.json": "membrane thickness",
    "TotalDensity.json": "total density profile",
    "FormFactor.json": "X-ray scattering form factor",
    "eq_times.json": "equilibration time",
}
TRAJECTORY_FORMATS = {
    ".xtc": "application/x-xtc",
    ".trr": "application/x-trr",
    ".dcd": "application/x-dcd",
    ".nc": "application/x-netcdf",
}

# Water is in every system and says nothing about it; ions do carry meaning.
KEYWORD_SKIP = frozenset({"SOL"})

# The generated block carries no marker comment: it is found by its own
# top-level key, always the last thing in the file. LEGACY_SENTINEL is kept only
# so an older banner is stripped along with the block rather than stranded.
BLOCK_KEY = "bioschema_properties:"
BLOCK_KEY_RE = re.compile(rf"^{re.escape(BLOCK_KEY)}", re.MULTILINE)
LEGACY_SENTINEL = "# --- bioschema_properties:"

BLOCK_ORDER = [
    "name", "alternateName", "description", "sameAs", "datePublished",
    "license", "articleLicense", "publisher", "version", "creator", "citation",
    "keywords", "measurementTechnique", "variableMeasured", "distribution",
    "isBasedOn", "isPartOf", "accessRights", "_source",
]

DOI_RE = re.compile(r"10\.\d{4,9}/[^\s;,\"']+")
ORCID_RE = re.compile(r"(\d{4}-\d{4}-\d{4}-\d{3}[\dX])", re.IGNORECASE)
DATE_RE = re.compile(r"^\d{4}(-\d{2}(-\d{2})?)?$")
NULLISH = frozenset({"null", "none", "na", "n/a", "??", "?", ""})
TRAILING = ".,;:)]}’\"'"


def _retry_delay(error, attempt):
    """Seconds to wait before the next attempt.

    Prefers a server-provided ``Retry-After`` header and otherwise falls back to
    exponential backoff with jitter, so concurrent callers do not retry in step.
    """
    headers = getattr(error, "headers", None)
    retry_after = headers.get("Retry-After") if headers is not None else None
    if retry_after:
        try:
            return min(float(retry_after), MAX_BACKOFF)
        except (TypeError, ValueError):
            pass
    return min(BACKOFF_BASE * (2**attempt), MAX_BACKOFF) + random.uniform(0, 0.5)


def fetch_json(url, timeout=DEFAULT_TIMEOUT):
    """GET a JSON document, or ``None`` when the resource is unavailable.

    Transient failures are retried; definitive ones (notably 404, meaning the
    DOI is not registered with this agency) are not.
    """
    request = urllib.request.Request(url)  # noqa: S310 - https URLs built above
    request.add_header("User-Agent", USER_AGENT)

    for attempt in range(MAX_RETRIES + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as error:
            if error.code not in RETRYABLE_STATUS or attempt == MAX_RETRIES:
                if error.code != 404:
                    print(f"  warning: {url} -> HTTP {error.code}")
                return None
            time.sleep(_retry_delay(error, attempt))
        except (urllib.error.URLError, TimeoutError, ValueError) as error:
            if attempt == MAX_RETRIES:
                print(f"  warning: {url} -> {error}")
                return None
            time.sleep(_retry_delay(error, attempt))
    return None


# ---------------------------------------------------------------------------
# SPDX licence resolution
# ---------------------------------------------------------------------------


def canonical_license_url(url):
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

    def resolve(self, uri):
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


def dataset_license(spdx):
    """This repository's own licence, resolved through the SPDX list."""
    if spdx is not None:
        licence, _ = spdx.resolve(DATASET_LICENSE_URI)
        if licence:
            return license_block(licence, DATASET_LICENSE_URI)
    return {"spdx": DATASET_LICENSE_SPDX, "url": DATASET_LICENSE_URI}


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


def access_rights(rights_list):
    """Keep the COAR access status ``resolve_license`` skips, rather than lose it."""
    for entry in rights_list or []:
        uri = (entry.get("rightsUri") or "").strip()
        if uri.startswith(ACCESS_RIGHTS_PREFIX):
            return uri[len(ACCESS_RIGHTS_PREFIX):]
    return None


# ---------------------------------------------------------------------------
# Field normalisation
# ---------------------------------------------------------------------------


def normalize_doi(value):
    """Extract a bare DOI from a URL, a ``doi:`` prefix or surrounding prose."""
    match = DOI_RE.search(str(value or ""))
    return match.group(0).rstrip(TRAILING) if match else None


DOI_ROLES = ("lookup", "cited", "article", "deposition")
Dois = namedtuple("Dois", DOI_ROLES)


def record_dois(readme, kind):
    """The record's DOIs, split by the role each one plays.

    One record can carry two DOIs that mean different things, and they are not
    interchangeable, so they are resolved once here rather than re-picked at
    each use site:

    ``lookup``      which registry record is fetched. For an experiment that is
                    the article whenever there is one: CrossRef carries the
                    authors, the journal and the publication date, and a raw
                    data deposition usually carries none of them.
    ``cited``       what the record points at as the thing to cite -- its
                    ``citation`` entry, its ``sameAs`` and the closing sentence
                    of its description. Data first, as
                    ``docs/src/schemas/experiment_metadata.md`` states: cite the
                    data, and fall back to the article only when it is the sole
                    DOI.
    ``article``     the parent work, used by ``part_of()``. An article is a
                    different relation from a citation: the dataset is *part of*
                    the paper it was digitised from while still *citing* the
                    deposition holding the raw values.
    ``deposition``  the raw data deposition, the parent work for a simulation
                    and for an experiment that has no article.

    ``DOI`` is the deprecated experiment spelling of ``ARTICLE_DOI``; 29 records
    still carry it, and it also holds ``unpublished/<slug>`` values, which
    ``normalize_doi`` rejects.
    """
    if kind == "simulations":
        # A simulation has one DOI, the Zenodo deposition holding the
        # trajectory. Nothing is cited from it: the deposition *is* this record,
        # so it belongs in isPartOf and distribution, not in a citation list.
        deposition = normalize_doi(readme.get("DOI"))
        return Dois(lookup=deposition, cited=deposition,
                    article=None, deposition=deposition)

    article = normalize_doi(readme.get("ARTICLE_DOI") or readme.get("DOI"))
    deposition = normalize_doi(readme.get("DATA_DOI"))
    return Dois(lookup=article or deposition, cited=deposition or article,
                article=article, deposition=deposition)


def normalize_orcid(value):
    match = ORCID_RE.search(str(value or ""))
    return f"https://orcid.org/{match.group(1).upper()}" if match else None


def iso_date(value):
    """Longest valid ISO prefix: ``YYYY-MM-DD``, ``YYYY-MM`` or ``YYYY``.

    Registry dates are usually complete but degrade to a year, and storing a
    year is better than storing nothing.
    """
    if value is None:
        return None
    text = str(value).strip()
    for pattern in (r"^\d{4}-\d{2}-\d{2}", r"^\d{4}-\d{2}", r"^\d{4}"):
        match = re.match(pattern, text)
        if match:
            return match.group(0)
    return None


def clean_text(value):
    """Collapse whitespace in an API-supplied string."""
    if value is None:
        return None
    return re.sub(r"\s+", " ", str(value)).strip() or None


def number(value):
    """Format a measured quantity without a spurious trailing ``.0``.

    Temperatures are written both as ``298`` and ``298.0`` across the corpus and
    the two mean the same thing; rendering them differently would give one system
    two titles.
    """
    try:
        return f"{float(value):g}"
    except (TypeError, ValueError):
        return clean_text(value)


def strip_markup(value):
    """Drop markup from an abstract.

    Zenodo descriptions contain HTML and CrossRef abstracts are JATS XML;
    neither belongs in a JSON-LD description, which is plain text.
    """
    if value is None:
        return None
    text = re.sub(r"<[^>]+>", " ", str(value))
    for entity, char in (("&amp;", "&"), ("&lt;", "<"), ("&gt;", ">"),
                         ("&quot;", '"'), ("&#39;", "'"), ("&nbsp;", " ")):
        text = text.replace(entity, char)
    return clean_text(text)


def parse_publication(value):
    """Turn a legacy ``PUBLICATION`` field into citation entries.

    These were hand-entered and inconsistent: bare DOIs, doi.org URLs, several
    DOIs joined by semicolons, and free-text references that sometimes carry a
    DOI in parentheses and sometimes none at all. ``schema.org/citation``
    accepts Text as well as CreativeWork, so a reference with no DOI is kept
    verbatim rather than dropped.
    """
    if value is None:
        return []
    text = str(value).strip()
    if text.lower() in NULLISH:
        return []
    entries = []
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        dois = DOI_RE.findall(part)
        if dois:
            entries.extend(d.rstrip(TRAILING) for d in dois)
        else:
            entries.append(part)
    return entries


def dedupe(entries):
    seen, unique = set(), []
    for entry in entries:
        if entry is None:
            continue
        key = str(entry).lower()
        if key not in seen:
            seen.add(key)
            unique.append(entry)
    return unique


# ---------------------------------------------------------------------------
# Registry mapping
# ---------------------------------------------------------------------------


def _publisher_name(publisher):
    """DataCite 4.5 allows publisher to be an object rather than a string."""
    return publisher.get("name") if isinstance(publisher, dict) else publisher


def source_block(api, doi, licence, spdx):
    """Provenance for the generated block.

    Records the SPDX list version whenever a licence actually resolved, so a
    surprising identifier can be traced back to the list that produced it --
    SPDX deprecates and re-scopes identifiers between releases.
    """
    source = {"api": api, "doi": doi}
    if spdx is not None and licence and licence.get("spdx"):
        source["spdxLicenseList"] = spdx.version
    return source


def from_datacite(payload, doi, readme, spdx):
    """Map a DataCite record onto Bioschemas properties."""
    attributes = (payload or {}).get("data", {}).get("attributes", {})
    if not attributes:
        return None, []

    notes = []
    titles = attributes.get("titles") or []
    dates = {d.get("dateType"): d.get("date") for d in attributes.get("dates") or []}
    published = iso_date(
        dates.get("Issued") or dates.get("Published") or dates.get("Available")
        or attributes.get("publicationYear")
    )

    licence, ambiguous = resolve_license(attributes.get("rightsList"), spdx)
    if ambiguous:
        notes.append("licence URI was ambiguous in the SPDX list")
    if licence is None:
        notes.append("no licence in DataCite rightsList")
    elif licence.get("spdx") is None:
        notes.append(f"licence did not resolve to SPDX: {licence.get('url')}")

    creators = []
    for creator in attributes.get("creators") or []:
        name = clean_text(creator.get("name")) or clean_text(
            " ".join(p for p in (creator.get("givenName"), creator.get("familyName")) if p)
        )
        if not name:
            continue
        entry = {"name": name}
        for identifier in creator.get("nameIdentifiers") or []:
            if (identifier.get("nameIdentifierScheme") or "").upper() == "ORCID":
                orcid = normalize_orcid(identifier.get("nameIdentifier"))
                if orcid:
                    entry["identifier"] = orcid
                    break
        creators.append(entry)
    if not creators:
        notes.append("no creators in DataCite record")

    # PUBLICATION is retired in favour of this list, so existing citations are
    # the authority. Re-deriving solely from PUBLICATION would erase them once
    # that field is gone.
    citations = [str(c) for c in ((readme.get("bioschema_properties") or {}).get("citation") or [])]
    citations += parse_publication(readme.get("PUBLICATION"))
    for related in attributes.get("relatedIdentifiers") or []:
        if (related.get("relationType") in CITATION_RELATIONS
                and (related.get("relatedIdentifierType") or "").upper() == "DOI"):
            citations.append(normalize_doi(related.get("relatedIdentifier")))

    block = {
        # Stashed for part_of(): this names the deposition, not the record, and
        # several hundred records share one deposition.
        "_article_title": strip_markup(titles[0].get("title")) if titles else None,
        "datePublished": published,
        "license": licence,
        "publisher": clean_text(_publisher_name(attributes.get("publisher"))),
        "version": clean_text(attributes.get("version")),
        "creator": creators,
        "citation": dedupe(citations),
    }
    rights = access_rights(attributes.get("rightsList"))
    if rights:
        block["accessRights"] = rights
    block["_subjects"] = [s.get("subject") for s in attributes.get("subjects") or [] if s.get("subject")]
    block["_source"] = source_block("datacite", doi, licence, spdx)
    return block, notes


def from_crossref(payload, doi, readme, spdx):
    """Map a CrossRef record onto Bioschemas properties."""
    message = (payload or {}).get("message") or {}
    if not message:
        return None, []

    notes = []
    titles = message.get("title") or []
    issued = (message.get("issued", {}).get("date-parts") or [[]])[0]
    published = iso_date("-".join(f"{p:02d}" if i else str(p) for i, p in enumerate(issued))) if issued else None

    licence, ambiguous = None, False
    entries = message.get("license") or []
    for entry in [e for e in entries if e.get("content-version") == "vor"] or entries:
        uri = (entry.get("URL") or "").strip()
        if not uri:
            continue
        if spdx is not None:
            resolved, ambiguous = spdx.resolve(uri)
            if resolved:
                licence = license_block(resolved, uri)
                break
        if licence is None:
            licence = {"spdx": None, "name": None, "url": uri}
    if licence is None:
        notes.append("no licence in CrossRef record")
    elif licence.get("spdx") is None:
        notes.append(f"licence did not resolve to SPDX: {licence.get('url')}")

    creators = []
    for author in message.get("author") or []:
        name = clean_text(" ".join(p for p in (author.get("given"), author.get("family")) if p)) or clean_text(
            author.get("name")
        )
        if not name:
            continue
        entry = {"name": name}
        orcid = normalize_orcid(author.get("ORCID"))
        if orcid:
            entry["identifier"] = orcid
        creators.append(entry)
    if not creators:
        notes.append("no authors in CrossRef record")

    existing = [str(c) for c in ((readme.get("bioschema_properties") or {}).get("citation") or [])]
    block = {
        # Stashed for part_of(): this names the article, not the record. It goes
        # through strip_markup because CrossRef titles carry JATS <sup>/<i>.
        "_article_title": strip_markup(titles[0]) if titles else None,
        "datePublished": published,
        "license": licence,
        "publisher": clean_text(message.get("publisher")),
        "version": None,
        "creator": creators,
        # The DOI this record cites is *not* necessarily the one fetched here:
        # which one it is follows record_dois(), and enrich() applies it.
        "citation": dedupe(existing + parse_publication(readme.get("PUBLICATION"))),
    }
    containers = message.get("container-title") or []
    if containers:
        # Stashed for enrich(): the journal is the parent of the *article*, not
        # of this dataset, so it nests inside isPartOf rather than replacing it.
        block["_journal"] = clean_text(containers[0])
    block["_source"] = source_block("crossref", doi, licence, spdx)
    return block, notes


def with_dataset_license(block, licence):
    """Move a fetched licence to ``articleLicense`` and set the dataset licence.

    Rebuilt rather than mutated so the two stay adjacent and ordered.
    """
    article = block.get("license")
    rebuilt = {}
    for key, value in block.items():
        if key == "license":
            rebuilt["license"] = licence
            if article:
                rebuilt["articleLicense"] = article
        else:
            rebuilt[key] = value
    rebuilt.setdefault("license", licence)
    return rebuilt


# ---------------------------------------------------------------------------
# Local enrichment
# ---------------------------------------------------------------------------


def edam(code, name):
    return {
        "@type": "DefinedTerm",
        "name": name,
        "termCode": code,
        "inDefinedTermSet": EDAM_SET,
        "url": f"{EDAM_SET}/{code}",
    }


def obo(code, name):
    return {
        "@type": "DefinedTerm",
        "name": name,
        "termCode": code,
        "inDefinedTermSet": f"{OBO_SET}/chmo.owl",
        "url": f"{OBO_SET}/{code}",
    }


def data_root_of(path):
    """Walk up from a README.yaml to the databank root."""
    for parent in Path(path).resolve().parents:
        if (parent / "Molecules").is_dir() and (parent / "Simulations").is_dir():
            return parent
    # Both layouts put README.yaml five levels below the root.
    return Path(path).resolve().parents[5]


def record_kind(path):
    return "simulations" if "Simulations" in Path(path).resolve().parts else "experiments"


def experiment_kind(path):
    return "xray" if "FormFactors" in Path(path).resolve().parts else "nmr"


def molecule_names(data_root):
    """Databank molecule id -> chemical name, from ``Molecules/*/*/metadata.yaml``."""
    names = {}
    for meta in Path(data_root).glob("Molecules/*/*/metadata.yaml"):
        try:
            block = yaml.safe_load(meta.read_text(encoding="utf-8")) or {}
        except (yaml.YAMLError, OSError):
            continue
        entry = block.get("NMRlipids") or {}
        if entry.get("name"):
            names[meta.parent.name] = entry["name"]
    return names


def split_subjects(subjects):
    """Zenodo subjects are free text; several pack a list into one string."""
    out = []
    for subject in subjects or []:
        for part in re.split(r"[,;]", str(subject)):
            part = part.strip()
            if part and len(part) < 80:
                out.append(part)
    return out


def dedupe_keywords(entries):
    seen, out = set(), []
    for entry in entries:
        key = (entry["name"] if isinstance(entry, dict) else str(entry)).lower()
        if key and key not in seen:
            seen.add(key)
            out.append(entry)
    return out


def _molecule_count(entry):
    """Total molecule count from a simulation ``COMPOSITION`` entry.

    ``COUNT`` is either a scalar or one entry per leaflet, and a leaflet entry is
    itself sometimes a list (united-atom records split a residue across lines).
    """
    count = (entry or {}).get("COUNT")
    if isinstance(count, list):
        total = 0
        for item in count:
            total += sum(item) if isinstance(item, list) else (item or 0)
        return float(total)
    return float(count or 0)


def composition_items(readme, kind):
    """Membrane components as ``(id, amount)``, largest share first.

    Amounts are molecule counts for simulations and molar fractions for
    experiments. ``MOLAR_FRACTIONS`` is the deprecated spelling of
    ``MEMBRANE_COMPOSITION`` and is the only composition 22 experiment records
    have, so it is read rather than leaving those records with no composition at
    all. Water is skipped: it is in every system and says nothing about it.
    """
    if kind == "simulations":
        raw = {mol: _molecule_count(entry)
               for mol, entry in (readme.get("COMPOSITION") or {}).items()}
    else:
        raw = readme.get("MEMBRANE_COMPOSITION") or readme.get("MOLAR_FRACTIONS") or {}

    items = []
    for mol, amount in raw.items():
        if mol in KEYWORD_SKIP:
            continue
        try:
            items.append((mol, float(amount)))
        except (TypeError, ValueError):
            items.append((mol, 0.0))
    # Descending share, then alphabetical, so a rerun cannot reorder the title.
    items.sort(key=lambda pair: (-pair[1], pair[0]))
    return items


def format_ratio(items):
    """``POPC`` for one component, ``POPC/POPE (95:5)`` for a mixture."""
    if not items:
        return "lipid"
    if len(items) == 1:
        return items[0][0]
    total = sum(amount for _, amount in items) or 1.0
    shares = ":".join(f"{100 * amount / total:.0f}" for _, amount in items)
    return "/".join(mol for mol, _ in items) + f" ({shares})"


def solution_ids(readme):
    """Ions actually present, by databank id. Zero-valued entries are padding."""
    out = []
    for mol, amount in (readme.get("SOLUTION_COMPOSITION") or {}).items():
        try:
            present = float(amount) != 0
        except (TypeError, ValueError):
            present = bool(amount)
        if present:
            out.append(mol)
    return sorted(out)


def composition_keywords(readme, kind, names):
    ids = [mol for mol, _ in composition_items(readme, kind)]
    if kind != "simulations":
        ids += solution_ids(readme)
    words = []
    for mol in ids:
        words.append(mol)
        if names.get(mol):
            words.append(names[mol])
    return words


def hydration_phrase(readme):
    """How wet the sample is, or ``None`` when the record does not say.

    ``TOTAL_HYDRATION`` is water mass %. Where it is missing the deprecated
    ``TOTAL_LIPID_CONCENTRATION`` is read the way ``Experiment.get_hydration``
    reads it -- as a lipid molarity converted to waters per lipid against water's
    own 55.5 M -- because for three records it is the only thing distinguishing
    the members of a hydration series.
    """
    hydration = readme.get("TOTAL_HYDRATION")
    if hydration is not None:
        return f"{number(hydration)}% water"

    lipid = readme.get("TOTAL_LIPID_CONCENTRATION")
    if lipid is None:
        return None
    if str(lipid).strip().lower() == "full hydration":
        return "full hydration"
    try:
        return f"{55.5 / float(lipid):.0f} waters per lipid"
    except (TypeError, ValueError, ZeroDivisionError):
        return None


def chmo_for_method(readme, path):
    """Map the record onto the most specific CHMO term.

    The directory decides NMR vs X-ray, not the presence of an ``NMR:`` block:
    around twenty OrderParameters records carry no NMR block, and keying off
    that would file them under small-angle X-ray scattering.
    """
    if experiment_kind(path) == "xray":
        return CHMO_SAXS
    method = str((readme.get("NMR") or {}).get("METHOD") or "")
    if method.startswith("PDLF:R"):
        return CHMO_PDLF_R
    if method.startswith("PDLF"):
        return CHMO_PDLF
    if method.startswith(("2H", "CDLF")):
        return CHMO_SSNMR
    return CHMO_NMR


def measurement_technique(readme, path, kind):
    """Ontology terms first, then the free-text detail no term captures."""
    if kind == "simulations":
        out = [edam(*EDAM_MD_OP)]
        engine, version = readme.get("SOFTWARE"), readme.get("SOFTWARE_VERSION")
        if engine:
            out.append(f"Molecular dynamics simulation ({engine}" + (f" {version}" if version else "") + ")")
        if readme.get("FF"):
            out.append(f"{readme['FF']} force field")
        return out

    nmr, xray = readme.get("NMR") or {}, readme.get("XRAY") or {}
    out = [obo(*chmo_for_method(readme, path))]
    if nmr or experiment_kind(path) == "nmr":
        out.append(edam(*EDAM_NMR_TOPIC))
        if nmr.get("METHOD"):
            out.append(f"NMR: {nmr['METHOD']}")
        if nmr.get("INSTRUMENT"):
            out.append(nmr["INSTRUMENT"])
    else:
        out.append(edam(*EDAM_XRAY_TOPIC))
        label = "X-ray scattering"
        if xray.get("SAMPLE_TYPE"):
            label += f" ({xray['SAMPLE_TYPE']})"
        out.append(label)
        if xray.get("SOURCE"):
            out.append(xray["SOURCE"])
    return out


def variables_measured(path, kind):
    """Read off which analyses actually exist beside this README."""
    directory = Path(path).resolve().parent
    found = [label for name, label in VARIABLES.items() if (directory / name).is_file()]
    if any(directory.glob("*OrderParameters.json")):
        found.append("C-H bond order parameter")
    if any(directory.glob("*_FormFactor.json")):
        found.append("X-ray scattering form factor")
    if kind == "experiments" and not found:
        found.append("C-H bond order parameter" if experiment_kind(path) == "nmr"
                     else "X-ray scattering form factor")

    seen, out = set(), []
    for item in found:
        if item in seen:
            continue
        seen.add(item)
        entry = {"@type": "PropertyValue", "name": item}
        code, symbol = UNITS.get(item, (None, None))
        if code:
            entry["unitCode"] = code
        if symbol:
            entry["unitText"] = symbol
        out.append(entry)
    return out


def distribution(readme, path, kind, doi):
    """Where the data actually lives."""
    if kind == "simulations":
        files = []
        for key in ("TRJ", "TPR", "TOP", "GRO", "PSF", "PDB"):
            for group in readme.get(key) or []:
                files += [f for f in (group or []) if isinstance(f, str)]
        entry = {"@type": "DataDownload"}
        if doi:
            entry["contentUrl"] = f"https://doi.org/{doi}"
        entry["name"] = readme.get("SYSTEM") or "trajectory deposition"
        if readme.get("TRAJECTORY_SIZE"):
            entry["contentSize"] = readme["TRAJECTORY_SIZE"]
        for candidate in files:
            suffix = Path(candidate).suffix.lower()
            if suffix in TRAJECTORY_FORMATS:
                entry["encodingFormat"] = TRAJECTORY_FORMATS[suffix]
                break
        if files:
            entry["hasPart"] = files
        return [entry]

    # Experiment data lives in this repository as JSON beside the README. No
    # contentUrl: the public URL is deployment-specific, so the web frontend
    # supplies it.
    directory = Path(path).resolve().parent
    return [
        {"@type": "DataDownload", "name": f.name, "encodingFormat": "application/json"}
        for f in sorted(directory.glob("*.json"))
    ]


def is_based_on(readme, kind):
    """Simulations record which experiment datasets they are compared against."""
    if kind != "simulations":
        return []
    experiment = readme.get("EXPERIMENT") or {}
    refs = [f"experiments/OrderParameters/{e}"
            for entries in (experiment.get("ORDERPARAMETER") or {}).values()
            for e in entries or []]
    refs += [f"experiments/FormFactors/{e}" for e in experiment.get("FORMFACTOR") or []]
    return dedupe(refs)


def part_of(block, dois, kind):
    """What this dataset is part of.

    For an experiment, preferably the article the values were digitised from,
    with the journal nested one level down as the article's own parent. Where
    there is no article -- the nmrXiv records -- the deposited dataset serves
    instead. For a simulation it is the Zenodo deposition holding the trajectory.
    Records with neither get nothing rather than an invented parent.

    The article wins here even where a ``DATA_DOI`` outranks it as the DOI to
    cite, and the two rules do not conflict: being *part of* a paper and citing
    the deposition that holds the raw values are different relations, and a
    record carrying both DOIs states both -- ``isPartOf`` the article,
    ``citation`` the deposition.

    This is also where the registry-supplied title lands. It names the parent
    work, which is what it was always describing; the record's own ``name`` is
    composed from the record's own values.
    """
    source_doi = (block.get("_source") or {}).get("doi")
    title = block.get("_article_title")

    if kind == "simulations":
        # The Zenodo deposition is a real parent: it holds the trajectory, and up
        # to 27 records share one. It is also where the fetched title belongs,
        # now that the record's own name is composed rather than borrowed.
        deposition = dois.deposition
        if not deposition:
            return None
        entry = {
            "@type": "Dataset",
            "@id": f"https://doi.org/{deposition}",
            "identifier": deposition,
            "url": f"https://doi.org/{deposition}",
        }
        if source_doi == deposition:
            if title:
                entry["name"] = title
            if block.get("publisher"):
                entry["publisher"] = block["publisher"]
        return entry

    article = dois.article
    if article:
        entry = {
            "@type": "ScholarlyArticle",
            "@id": f"https://doi.org/{article}",
            "identifier": article,
            "url": f"https://doi.org/{article}",
        }
        if title and source_doi == article:
            entry["name"] = title
        if block.get("_journal"):
            entry["isPartOf"] = {"@type": "Periodical", "name": block["_journal"]}
        return entry

    deposition = dois.deposition
    if deposition:
        entry = {
            "@type": "Dataset",
            "@id": f"https://doi.org/{deposition}",
            "identifier": deposition,
            "url": f"https://doi.org/{deposition}",
        }
        if source_doi == deposition:
            if title:
                entry["name"] = title
            if block.get("publisher"):
                entry["publisher"] = block["publisher"]
        return entry
    return None


def measured_quantity(path):
    """What an experiment record actually holds."""
    return ("X-ray scattering form factor" if experiment_kind(path) == "xray"
            else "C-H bond order parameters")


def short_technique(readme, path):
    """A title-length technique label, from the same fields as ``chmo_for_method``."""
    if experiment_kind(path) == "xray":
        sample = (readme.get("XRAY") or {}).get("SAMPLE_TYPE")
        return f"SAXS, {sample}" if sample else "SAXS"
    method = str((readme.get("NMR") or {}).get("METHOD") or "")
    if method.startswith("PDLF:R"):
        return "R-PDLF NMR"
    if method.startswith("PDLF"):
        return "PDLF NMR"
    if method.startswith("2H"):
        return "2H NMR"
    if method.startswith("CDLF"):
        return "CDLF NMR"
    return "NMR"


def surname(name):
    """Family name from a creator entry.

    DataCite writes ``Family, Given`` and CrossRef ``Given Family``, so the comma
    decides which end to take. Compound names are kept whole on the CrossRef side
    only where the parts are lowercase particles (``van der Berg``).
    """
    text = clean_text(name)
    if not text:
        return None
    if "," in text:
        return text.split(",")[0].strip() or None
    parts = text.split()
    if not parts:
        return None
    start = len(parts) - 1
    while start > 0 and parts[start - 1][:1].islower():
        start -= 1
    return " ".join(parts[start:])


def source_tag(block, readme, path, kind):
    """The trailing bracket that keeps two similar records apart.

    Simulations use their databank ``ID``, which BilayerData's own ``CheckIDs.sh``
    keeps unique. Experiments have no such field, so they use the first author and
    year -- which distinguishes the systems that two different groups measured
    independently -- falling back to the record's own directory when there is no
    publication to name.
    """
    if kind == "simulations":
        identifier = readme.get("ID")
        return f"NMRlipids simulation {identifier}" if identifier is not None else None

    creators = block.get("creator") or []
    family = surname(creators[0].get("name")) if creators else None
    year = iso_date(block.get("datePublished"))
    if family:
        return f"{family} {year[:4]}" if year else family

    parts = Path(path).resolve().parts
    # .../experiments/<kind>/unpublished/<slug>/<index>/README.yaml
    if "unpublished" in parts:
        return f"unpublished, {parts[-3]}"
    return None


def compose_name(readme, path, kind, block):
    """The record's title, pasted together from the record's own values."""
    items = composition_items(readme, kind)
    system = f"{format_ratio(items)} bilayer"
    temperature = number(readme.get("TEMPERATURE"))

    if kind == "simulations":
        title = f"Molecular dynamics trajectory of a {system}"
        if temperature:
            title += f" at {temperature} K"
        detail = []
        if readme.get("FF"):
            detail.append(clean_text(readme["FF"]))
        engine = clean_text(readme.get("SOFTWARE"))
        if engine:
            version = clean_text(readme.get("SOFTWARE_VERSION"))
            detail.append(engine.upper() + (f" {version}" if version else ""))
        if readme.get("TRJLENGTH"):
            detail.append(f"{float(readme['TRJLENGTH']) / 1000:.0f} ns")
        if detail:
            title += " (" + ", ".join(detail) + ")"
    else:
        title = f"{measured_quantity(path)} of a {system}"
        if temperature:
            title += f" at {temperature} K"
        hydration = hydration_phrase(readme)
        if hydration:
            title += f", {hydration}"
        ions = solution_ids(readme)
        additives = sorted((readme.get("ADDITIONAL_MOLECULES") or {}).keys())
        if ions:
            title += ", with " + "/".join(ions)
        if additives:
            title += (" and " if ions else ", with ") + "/".join(additives)
        title += f" ({short_technique(readme, path)})"

    tag = source_tag(block, readme, path, kind)
    return title + (f" [{tag}]" if tag else "")


def _named(mol, names):
    return f"{mol} ({names[mol]})" if names.get(mol) else mol


def _component(mol, names, share=None, count=None):
    """One component of a bilayer: ``200 POPC (1-palmitoyl-..., 90 mol%)``.

    The chemical name and the share share one parenthesis; keeping them apart
    reads as two unrelated asides.
    """
    inner = [names[mol]] if names.get(mol) else []
    if share is not None:
        inner.append(f"{share:.0f} mol%")
    phrase = f"{count:.0f} {mol}" if count is not None else mol
    return phrase + (" (" + ", ".join(inner) + ")" if inner else "")


def _listed(entries):
    if len(entries) == 1:
        return entries[0]
    return ", ".join(entries[:-1]) + " and " + entries[-1]


def compose_description(readme, path, kind, block, dois, names):
    """Three sentences: what the system is, how it was measured, where it came from."""
    items = composition_items(readme, kind)
    temperature = number(readme.get("TEMPERATURE"))

    if kind == "simulations":
        total = sum(amount for _, amount in items) or 1.0
        parts = [_component(mol, names, count=amount,
                            share=None if len(items) == 1 else 100 * amount / total)
                 for mol, amount in items]
        system = ("Molecular dynamics trajectory of a lipid bilayer of "
                  + (_listed(parts) if parts else "unrecorded composition"))
        if temperature:
            system += f" at {temperature} K"

        force_field = clean_text(readme.get("FF"))
        engine = clean_text(readme.get("SOFTWARE"))
        if engine:
            version = clean_text(readme.get("SOFTWARE_VERSION"))
            engine = engine.upper() + (f" {version}" if version else "")
        second = "Simulated"
        if force_field:
            second += f" with {force_field}"
        if engine:
            second += f" in {engine}"
        if not force_field and not engine:
            # No record says what produced it, so do not claim a method.
            second = "Trajectory run"
        if readme.get("TRJLENGTH"):
            second += f" for {float(readme['TRJLENGTH']) / 1000:.0f} ns"
        if readme.get("NUMBER_OF_ATOMS"):
            second += f" ({readme['NUMBER_OF_ATOMS']} atoms)"

        identifier = readme.get("ID")
        doi = dois.cited
        third = "Part of the NMRlipids Databank"
        if identifier is not None:
            third = f"Deposited as NMRlipids Databank simulation {identifier}"
        if doi:
            third += f" and available from https://doi.org/{doi}"
        return " ".join(f"{s}." for s in (system, second, third))

    total = sum(amount for _, amount in items) or 1.0
    parts = [_component(mol, names,
                        share=None if len(items) == 1 else 100 * amount / total)
             for mol, amount in items]
    system = (f"Experimental {measured_quantity(path)} for a lipid bilayer of "
              + (_listed(parts) if parts else "unrecorded composition"))
    if temperature:
        system += f" at {temperature} K"
    hydration = hydration_phrase(readme)
    if hydration == "full hydration":
        system += ", fully hydrated"
    elif hydration:
        system += f", hydrated to {hydration}"

    conditions = []
    ions = solution_ids(readme)
    if ions:
        conditions.append("ions " + _listed([_named(i, names) for i in ions]))
    additives = sorted((readme.get("ADDITIONAL_MOLECULES") or {}).keys())
    if additives:
        conditions.append("additives " + _listed(additives))
    ph = readme.get("PH")
    if ph is not None and str(ph).strip().lower() not in NULLISH | {"unknown"}:
        how = clean_text(readme.get("PH_METHOD"))
        conditions.append(f"pH {number(ph)}"
                          + (f" ({how})" if how and how.lower() != "unknown" else ""))
    if conditions:
        system += ", with " + "; ".join(conditions)

    second = f"Measured by {chmo_for_method(readme, path)[1]}"
    if experiment_kind(path) == "xray":
        xray = readme.get("XRAY") or {}
        if xray.get("SAMPLE_TYPE"):
            second += f" on {clean_text(xray['SAMPLE_TYPE'])} samples"
        if xray.get("SOURCE"):
            second += f" at {clean_text(xray['SOURCE'])}"
    else:
        nmr = readme.get("NMR") or {}
        if nmr.get("METHOD"):
            second += f" ({clean_text(nmr['METHOD'])})"
        if nmr.get("INSTRUMENT"):
            second += f" on a {clean_text(nmr['INSTRUMENT'])}"

    # The cited DOI, so the sentence, sameAs and citation all name one source.
    doi = dois.cited
    third = (f"Values digitised into the NMRlipids Databank from https://doi.org/{doi}"
             if doi else "Unpublished data contributed to the NMRlipids Databank")
    return " ".join(f"{s}." for s in (system, second, third))


def experiment_citations(block, dois):
    """Apply the data-first citation rule to one experiment block.

    ``docs/src/schemas/experiment_metadata.md``: a ``DATA_DOI`` is what gets
    cited, and the ``ARTICLE_DOI`` only when it is the sole DOI. Where both are
    given the article is actively removed rather than merely not added, so a
    block written before the rule existed -- every enriched experiment in
    BilayerData carries the article there -- is corrected on the next run. The
    article is not lost: ``part_of()`` states it as the parent work.

    Everything else already in the list stays. Those entries come from the
    legacy ``PUBLICATION`` field and from DataCite ``relatedIdentifiers``; they
    are related publications, which the rule says nothing about.
    """
    cites = [str(c) for c in (block.get("citation") or [])]
    if dois.deposition and dois.article:
        cites = [c for c in cites if normalize_doi(c) != dois.article]
    return dedupe(cites + [dois.cited])


def enrich(block, readme, path, kind, dois, names):
    """Add every Bioschemas property derivable from local data."""
    subjects = block.pop("_subjects", [])
    # Composed, never fetched: a registry title names the paper or deposition, and
    # hundreds of records share one of those.
    block["name"] = compose_name(readme, path, kind, block)
    block["description"] = compose_description(readme, path, kind, block, dois, names)
    if dois.cited:
        # `identifier` is deliberately not written: the web frontend derives it
        # from the record's own DOI field. sameAs keeps the resolvable link, and
        # it names the cited DOI rather than the one that was looked up: what
        # this dataset is a copy of is the data, not the paper describing it.
        block["sameAs"] = f"https://doi.org/{dois.cited}"
    if kind == "simulations" and readme.get("SYSTEM"):
        block["alternateName"] = readme["SYSTEM"]

    terms = EDAM_SIMULATION if kind == "simulations" else (
        EDAM_XRAY if experiment_kind(path) == "xray" else EDAM_NMR
    )
    block["keywords"] = dedupe_keywords(
        [edam(code, label) for code, label in terms]
        + composition_keywords(readme, kind, names)
        + split_subjects(subjects)
    )
    block["measurementTechnique"] = measurement_technique(readme, path, kind)
    block["variableMeasured"] = variables_measured(path, kind)
    if kind == "experiments":
        block["citation"] = experiment_citations(block, dois)

    dist = distribution(readme, path, kind, dois.lookup)
    if dist:
        block["distribution"] = dist
    based = is_based_on(readme, kind)
    if based:
        block["isBasedOn"] = based

    parent = part_of(block, dois, kind)
    block.pop("_journal", None)
    block.pop("_article_title", None)
    if parent:
        block["isPartOf"] = parent
    else:
        block.pop("isPartOf", None)
    return block


# ---------------------------------------------------------------------------
# Reading and writing README.yaml
# ---------------------------------------------------------------------------


def ordered(block):
    """Canonical key order; anything unlisted keeps its place at the end."""
    out = {k: block[k] for k in BLOCK_ORDER if k in block}
    out.update({k: v for k, v in block.items() if k not in out})
    return out


def prune(block):
    """Drop absent values, as ``Molecules/*/metadata.yaml`` does.

    An omitted property and one set to null mean the same thing in JSON-LD.
    Licences are the exception: an unresolved one keeps an explicit ``spdx:
    null`` beside the URI that failed to resolve, so the gap stays visible.
    """
    pruned = {}
    for key, value in block.items():
        if key in ("license", "articleLicense") and isinstance(value, dict):
            licence = {k: v for k, v in value.items() if v is not None}
            licence.setdefault("spdx", None)
            pruned[key] = licence
        elif value not in (None, [], {}):
            pruned[key] = value
    return pruned


def strip_existing_block(text):
    """Everything above the generated block.

    With no block present the file is returned byte for byte, gaining only the
    newline separating it from what gets appended: several experiment READMEs
    end in blank lines after the deprecated-field section, and normalising those
    away would be an edit to hand-written content.
    """
    index = text.find(LEGACY_SENTINEL)
    if index == -1:
        match = BLOCK_KEY_RE.search(text)
        index = match.start() if match else -1
    if index == -1:
        return text if text.endswith("\n") else text + "\n"
    return text[:index]


def render_block(block):
    return yaml.dump(
        {"bioschema_properties": block},
        sort_keys=False,
        allow_unicode=True,
        default_flow_style=False,
        width=4096,
    )


def volatile_stripped(block):
    """A block with ``retrieved`` removed, for change detection.

    Without this the retrieval date alone would make every rerun a diff.
    """
    if not isinstance(block, dict):
        return block
    copy = json.loads(json.dumps(block, sort_keys=True, default=str))
    if isinstance(copy.get("_source"), dict):
        copy["_source"].pop("retrieved", None)
    return copy


# Date fields normalised per record kind, mapped to whether a ``DD/MM/YYYY``
# value may be reinterpreted as ISO.
#
# Only DATEOFRUNNING may: add_simulation.py writes it with strftime("%d/%m/%Y"),
# so the convention is known. FF_DATE is hand-entered and the corpus proves both
# conventions are in use -- 14/10/2025 and 26/01/2021 can only be DD/MM, while
# 3/22/21 can only be MM/DD -- alongside values like "?/?/2020" and "10/2023".
# Reinterpreting those would silently corrupt them, so slash-formatted values
# are left exactly as they are; both fields still get quoted when already ISO,
# which is what the schemas need (unquoted, YAML yields a date object).
#
# A simulation's DATE is not declared in readme_yaml_schema.json at all, so it
# is left as YAML parses it rather than rewritten for tidiness.
DATE_FIELDS = {
    "simulations": {"DATEOFRUNNING": True, "FF_DATE": False},
    "experiments": {"DATE": False},
}


def normalize_dates(text, kind):
    """Rewrite the record's date fields as quoted ISO strings.

    Quoting matters: unquoted, YAML parses ``2011-07-26`` into a date object and
    the schemas require these fields to be strings. ``DATEOFRUNNING`` also
    arrives as ``DD/MM/YYYY`` from ``add_simulation.py``.
    """
    fields = "|".join(DATE_FIELDS[kind])
    lines = text.splitlines(keepends=True)
    changed = False
    for i, line in enumerate(lines):
        match = re.match(rf"^({fields}):[ \t]*['\"]?([^'\"#\n]+?)['\"]?[ \t]*(#.*)?$",
                         line.rstrip("\n"))
        if not match:
            continue
        key, raw, comment = match.group(1), match.group(2).strip(), match.group(3)
        slashed = re.fullmatch(r"(\d{2})/(\d{2})/(\d{4})", raw)
        if slashed and DATE_FIELDS[kind][key]:
            value = f"{slashed.group(3)}-{slashed.group(2)}-{slashed.group(1)}"
        else:
            iso = re.match(r"^(\d{4}-\d{2}-\d{2})", raw)
            if not iso:
                continue
            value = iso.group(1)
        new = f"{key}: '{value}'" + (f"  {comment}" if comment else "")
        new += "\n" if line.endswith("\n") else ""
        if new != line:
            lines[i] = new
            changed = True
    return "".join(lines), changed


def drop_publication(text, block):
    """Retire ``PUBLICATION`` once ``citation`` demonstrably covers it.

    Refuses to remove anything the citation list does not already carry, so a
    hand-written reference cannot be lost to a parsing slip.
    """
    doc = yaml.safe_load(text) or {}
    if "PUBLICATION" not in doc:
        return text, False
    cites = [str(c) for c in (block.get("citation") or [])]
    for item in parse_publication(doc.get("PUBLICATION")):
        if item not in cites:
            print(f"  keeping PUBLICATION: {item!r} is not in citation")
            return text, False
    kept = [ln for ln in text.splitlines(keepends=True) if not ln.startswith("PUBLICATION:")]
    return "".join(kept), True


def resolve_doi(doi, kind, cache_dir):
    """Fetch one DOI, memoised on disk. Returns ``(payload, api)``.

    Simulations are Zenodo depositions, so DataCite answers. Experiments are
    usually journal articles, so CrossRef answers -- but a few cite an nmrXiv or
    Zenodo deposition instead, and those are registered with DataCite, so fall
    back to it rather than skipping them.
    """
    slug = urllib.parse.quote(doi, safe="")
    cached = cache_dir / f"{slug}.json"
    if cached.is_file():
        payload, api = json.loads(cached.read_text(encoding="utf-8"))
        return payload, api

    quoted = urllib.parse.quote(doi, safe="/")
    order = ([("datacite", DATACITE_URL)] if kind == "simulations"
             else [("crossref", CROSSREF_URL), ("datacite", DATACITE_URL)])
    payload, api = None, order[0][0]
    for candidate, template in order:
        payload = fetch_json(template.format(doi=quoted))
        if payload is not None:
            api = candidate
            break

    cache_dir.mkdir(parents=True, exist_ok=True)
    cached.write_text(json.dumps([payload, api]), encoding="utf-8")
    return payload, api


def process(path, spdx, names, cache_dir, dry_run=False):
    """Enrich one README.yaml in place. Returns True when the file changed."""
    path = Path(path)
    original = path.read_text(encoding="utf-8")
    readme = yaml.safe_load(original) or {}
    kind = record_kind(path)

    dois = record_dois(readme, kind)
    doi = dois.lookup

    block, notes = None, []
    if doi:
        payload, api = resolve_doi(doi, kind, cache_dir)
        if payload is None:
            print(f"  warning: DOI not found: {doi}")
        else:
            mapper = from_datacite if api == "datacite" else from_crossref
            block, notes = mapper(payload, doi, readme, spdx)
    if block is None:
        # unpublished/<slug> records have nothing to look up, but they are still
        # data this repository distributes, so they carry the dataset licence.
        block = {"_source": {"api": None, "doi": None}}
    for note in notes:
        print(f"  note: {note}")

    block = enrich(block, readme, path, kind, dois, names)
    if kind == "experiments":
        block = with_dataset_license(block, dataset_license(spdx))
    block = ordered(prune(block))

    text, dates_changed = normalize_dates(original, kind)
    if kind == "simulations":
        text, _ = drop_publication(text, block)

    existing = readme.get("bioschema_properties")
    if volatile_stripped(existing) == volatile_stripped(block):
        block.setdefault("_source", {})["retrieved"] = (
            ((existing or {}).get("_source") or {}).get("retrieved") or date.today().isoformat()
        )
    else:
        block.setdefault("_source", {})["retrieved"] = date.today().isoformat()

    updated = strip_existing_block(text) + render_block(block)
    if updated == original:
        return False
    if not dry_run:
        path.write_text(updated, encoding="utf-8")
    print(f"{'Would update' if dry_run else 'Updated'} bioschema metadata in {path}"
          + (" (dates normalised)" if dates_changed else ""))
    return True


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def check(path, spdx, strict=False):
    """Validate a generated block. Returns ``(level, message)`` pairs.

    ERROR means the block is unusable; WARNING means it is well-formed but
    incomplete. Two things warn rather than fail: a licence outside SPDX is the
    expected outcome for experiments, whose CrossRef licences are publisher
    text-mining terms, and a missing block is normal for a freshly uploaded
    record the workflow has not reached yet.
    """
    readme = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    block = readme.get("bioschema_properties")
    level = "ERROR" if strict else "WARNING"
    if block is None:
        return [(level, f"{level}: {path}: missing bioschema_properties block")]
    if not isinstance(block, dict):
        return [("ERROR", f"ERROR: {path}: bioschema_properties is not a mapping")]

    found = []

    def error(message):
        found.append(("ERROR", f"ERROR: {path}: {message}"))

    def warn(message):
        found.append(("WARNING", f"WARNING: {path}: {message}"))

    resolved = (block.get("_source") or {}).get("doi") is not None
    # name and description are composed from the record itself, so unlike the
    # registry-supplied properties they are expected even without a resolved DOI.
    if not clean_text(block.get("name")):
        error("name is empty")
    if not clean_text(block.get("description")):
        error("description is empty")

    published = block.get("datePublished")
    if resolved and published is None:
        error("datePublished is missing")
    elif published is not None and not DATE_RE.match(str(published)):
        error(f"datePublished {published!r} is not YYYY[-MM[-DD]]")

    licence = block.get("license")
    if not isinstance(licence, dict):
        warn("no licence: the DOI record declares none")
    elif licence.get("spdx") is None:
        warn(f"licence did not resolve to SPDX (url {licence.get('url')!r})")
    elif spdx is not None and licence["spdx"] not in spdx.by_id:
        error(f"license.spdx {licence['spdx']!r} is not in the SPDX list")

    article = block.get("articleLicense")
    if isinstance(article, dict) and article.get("spdx") is None:
        warn(f"articleLicense is outside SPDX -- publisher terms ({article.get('url')!r})")

    creators = block.get("creator")
    if resolved and not creators:
        error("creator is empty")
    for index, creator in enumerate(creators or []):
        if not isinstance(creator, dict) or not clean_text(creator.get("name")):
            error(f"creator[{index}] has no name")
        elif creator.get("identifier") and not ORCID_RE.search(str(creator["identifier"])):
            error(f"creator[{index}].identifier {creator['identifier']!r} is not an ORCID URI")

    # The data-first rule, checked against the record's own DOI fields rather
    # than re-derived: a block written before the rule, or hand-edited since,
    # still cites the article, and --check is how a whole databank is swept for
    # the records a rerun has to visit.
    kind = record_kind(path)
    dois = record_dois(readme, kind)
    cites = {normalize_doi(c) for c in (block.get("citation") or [])}
    if kind == "experiments":
        if dois.cited and dois.cited not in cites:
            error(f"citation does not include {dois.cited}, which is the DOI this record cites")
        if dois.deposition and dois.article and dois.article in cites:
            error(f"citation carries the article {dois.article}; DATA_DOI {dois.deposition} outranks it")
    if dois.cited and normalize_doi(block.get("sameAs")) != dois.cited:
        error(f"sameAs {block.get('sameAs')!r} does not point at the cited DOI {dois.cited}")

    return found


def duplicate_names(paths):
    """Records sharing a ``name``. Returns ``(level, message)`` pairs.

    Titles are composed rather than fetched precisely so that no two records
    share one, and a catalogue listing them is unusable if two do. Only a run
    over the whole databank can prove that, so this is worth pointing at
    everything rather than at a pull request's changed files.
    """
    claimed = {}
    for path in paths:
        try:
            readme = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
        except (yaml.YAMLError, OSError):
            continue
        block = readme.get("bioschema_properties")
        name = clean_text(block.get("name")) if isinstance(block, dict) else None
        if name:
            claimed.setdefault(name, []).append(path)

    found = []
    for name, holders in sorted(claimed.items()):
        if len(holders) > 1:
            listed = ", ".join(str(p) for p in holders)
            found.append(("ERROR", f"ERROR: duplicate name {name!r} in {len(holders)} "
                                   f"records: {listed}"))
    return found


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", help="README.yaml paths to enrich")
    parser.add_argument("--cache", type=Path, default=None,
                        help="directory for cached API responses")
    parser.add_argument("--dry-run", action="store_true",
                        help="report what would change without writing")
    parser.add_argument("--check", action="store_true",
                        help="validate existing blocks instead of writing")
    parser.add_argument("--strict", action="store_true",
                        help="with --check, fail on a missing block rather than warn")
    args = parser.parse_args()

    paths = [Path(f) for f in args.files]
    missing = [p for p in paths if not p.is_file()]
    if missing:
        sys.exit(f"error: no such file: {missing[0]}")

    root = data_root_of(paths[0])
    cache_dir = args.cache or Path(os.environ.get("RUNNER_TEMP", root)) / ".cache" / "bioschema"
    spdx = load_spdx(cache_dir)
    if spdx is None:
        print("warning: SPDX licence list unavailable; licences will not resolve")

    if args.check:
        failed = set()
        for path in paths:
            for level, message in check(path, spdx, args.strict):
                print(message)
                if level == "ERROR":
                    failed.add(path)
        duplicates = duplicate_names(paths)
        for _, message in duplicates:
            print(message)
        print(f"\n{len(paths) - len(failed)} of {len(paths)} records valid"
              + (f", {len(failed)} failing" if failed else "")
              + (f", {len(duplicates)} duplicated names" if duplicates else ""))
        sys.exit(1 if failed or duplicates else 0)

    names = molecule_names(root)
    changed = sum(process(p, spdx, names, cache_dir, args.dry_run) for p in paths)
    print(f"\n{changed} of {len(paths)} records {'would be ' if args.dry_run else ''}updated")


if __name__ == "__main__":
    main()
