"""
Validating what a run wrote.

``check`` re-reads one record and reports what is missing, malformed or
inconsistent with the record's own fields; ``duplicate_names`` is the one test
that cannot be made record by record, since two records only clash when both are
looked at. Both are what ``--check`` reports, and neither modifies anything.
"""

from pathlib import Path

import yaml

from .constants import DATE_RE, ORCID_RE
from .fields import deprecated_keys, record_dois, record_kind
from .helpers import clean_text, normalize_doi


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

    kind = record_kind(path)

    # A deprecated key warns rather than fails: the block generated from it is
    # correct, and what is wrong is the record around it, which
    # experiment_schema.json will not accept until the key is renamed. --check is
    # therefore also the sweep that lists which records are still to be migrated.
    for key, replacement in deprecated_keys(readme, kind).items():
        warn(f"{key} is deprecated and not accepted by experiment_schema.json; "
             f"rename it to {replacement}")

    # The data-first rule, checked against the record's own DOI fields rather
    # than re-derived: a block written before the rule, or hand-edited since,
    # still cites the article, and --check is how a whole databank is swept for
    # the records a rerun has to visit.
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

    An experiment's title ends in its source tag precisely so that no two
    experiments share one, so a clash there is an ERROR. A simulation's title
    carries no tag -- its databank ``ID`` is assigned only after merge -- and
    similar simulations can share one; the website's JSON-LD identifies them by
    their page, so a clash there is only a WARNING. Only a run over the whole
    databank can find either, so this is worth pointing at everything rather
    than at a pull request's changed files.
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
            simulations = all(record_kind(p) == "simulations" for p in holders)
            level = "WARNING" if simulations else "ERROR"
            found.append((level, f"{level}: duplicate name {name!r} in {len(holders)} "
                                 f"records: {listed}"))
    return found
