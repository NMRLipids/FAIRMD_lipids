"""
Turning a DataCite or CrossRef answer into a Bioschemas block.

One function per registry, because the two describe a work differently: DataCite
answers for the Zenodo depositions every simulation points at, CrossRef for the
journal articles most experiments point at. Both are mapped onto the same set of
properties here, so the rest of the pipeline never sees which registry replied.

``resolve_doi`` is the one entry point that touches the network, and it caches
each answer on disk so a re-run over the databank asks once per DOI.
"""

import json
import urllib.parse

from .constants import (
    CITATION_RELATIONS,
    CROSSREF_URL,
    DATACITE_URL,
)
from .helpers import (
    clean_text,
    dedupe,
    fetch_json,
    iso_date,
    normalize_doi,
    normalize_orcid,
    parse_publication,
    strip_markup,
)
from .licenses import access_rights, compact_license, license_block, resolve_license


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


def with_dataset_license(block, licence, article=None):
    """State the databank's licence, and move the fetched one where it belongs.

    An experiment record digitises values out of a source: the source carries
    its own terms, while what this repository distributes is the digitised
    values, under its own licence. That licence is the data repository's to
    state, in the ``experiments/bioschema_properties.yaml`` template every
    experiment inherits ``license`` from, and it is passed in from there; the
    compaction in ``expsim_metadata.records`` then drops it from the record
    again, since it is inherited anyway. Where the fetched licence goes depends on
    what answered for it, not on which registry did: it is an ``articleLicense``
    only when the DOI that was looked up is the article the values came from.
    The four nmrXiv records give a ``DATA_DOI`` and no article, and their
    deposition's licence is a property of that deposition -- so it lands on
    ``isPartOf``, which is the deposition, rather than being labelled as the
    licence of an article that does not exist.

    Rebuilt rather than mutated so the two stay adjacent and ordered.
    """
    fetched = block.get("license")
    source_doi = (block.get("_source") or {}).get("doi")
    describes_article = bool(article) and source_doi == article

    rebuilt = {}
    for key, value in block.items():
        if key == "license":
            rebuilt["license"] = licence
            if fetched and describes_article:
                rebuilt["articleLicense"] = fetched
        else:
            rebuilt[key] = value
    rebuilt.setdefault("license", licence)

    parent = rebuilt.get("isPartOf")
    if (fetched and not describes_article and isinstance(parent, dict)
            and parent.get("identifier") == source_doi):
        nested = compact_license(fetched) if isinstance(fetched, dict) else fetched
        rebuilt["isPartOf"] = {**parent, "license": nested}
    return rebuilt


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

    # Only an answer is cached. A failed lookup is a fact about today -- an
    # outage, an exhausted retry budget -- and caching it would suppress every
    # later attempt until someone deleted the file by hand.
    if payload is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        cached.write_text(json.dumps([payload, api]), encoding="utf-8")
    return payload, api
