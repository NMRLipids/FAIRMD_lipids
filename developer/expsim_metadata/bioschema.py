"""
Building the Bioschemas ``Dataset`` block of one record.

Everything the block states about a record that is not simply copied from a
registry is decided here: which ontology term describes the record's own method,
which keywords the composition contributes, what the analysis outputs beside the
README amount to as a distribution, and which work the record is part of.

What the record's *kind* alone decides -- the EDAM keywords and technique terms,
an experiment's measured variable, the databank's licence -- is not written
here at all: the ``bioschema_properties.yaml`` templates in the data repository
state it once for every record below them (see ``expsim_metadata.inheritance``).

``enrich`` is the whole of it in order, and is what ``expsim_metadata.records`` calls
once the registry answer has been fetched.
"""

import re
from pathlib import Path

from .constants import (
    OBO_SET,
    TRAJECTORY_FORMATS,
    UNITS,
    VARIABLES,
)
from .descriptions import compose_description, compose_name
from .fields import chmo_for_method, composition_items, experiment_kind, solution_ids
from .helpers import dedupe, normalize_doi


def obo(code, name):
    return {
        "@type": "DefinedTerm",
        "name": name,
        "termCode": code,
        "inDefinedTermSet": f"{OBO_SET}/chmo.owl",
        "url": f"{OBO_SET}/{code}",
    }


def split_subjects(subjects):
    """Zenodo subjects are free text; several pack a list into one string."""
    out = []
    for subject in subjects or []:
        for part in re.split(r"[,;]", str(subject)):
            part = part.strip()
            if part and len(part) < 80:
                out.append(part)
    return out


def _keyword_key(entry):
    return ((entry.get("name") or "") if isinstance(entry, dict) else str(entry)).lower()


def dedupe_keywords(entries, inherited=()):
    """Drop repeated keywords, by name and regardless of case.

    Seeded with the keywords the record inherits, so a free-text Zenodo subject
    such as "Molecular dynamics" is not repeated beside the EDAM term of that
    name that every simulation already inherits.
    """
    seen = {_keyword_key(entry) for entry in inherited}
    out = []
    for entry in entries:
        key = _keyword_key(entry)
        if key and key not in seen:
            seen.add(key)
            out.append(entry)
    return out


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


def measurement_technique(readme, path, kind):
    """The record's own method: its CHMO term, then the free-text detail no term captures.

    The technique terms shared by the whole kind -- EDAM's molecular dynamics
    operation, its NMR and X-ray diffraction topics -- are inherited from the
    templates, not written per record.
    """
    if kind == "simulations":
        out = []
        engine, version = readme.get("SOFTWARE"), readme.get("SOFTWARE_VERSION")
        if engine:
            out.append(f"Molecular dynamics simulation ({engine}" + (f" {version}" if version else "") + ")")
        if readme.get("FF"):
            out.append(f"{readme['FF']} force field")
        return out

    nmr, xray = readme.get("NMR") or {}, readme.get("XRAY") or {}
    out = [obo(*chmo_for_method(readme, path))]
    if nmr or experiment_kind(path) == "nmr":
        if nmr.get("METHOD"):
            out.append(f"NMR: {nmr['METHOD']}")
        if nmr.get("INSTRUMENT"):
            out.append(nmr["INSTRUMENT"])
    else:
        label = "X-ray scattering"
        if xray.get("SAMPLE_TYPE"):
            label += f" ({xray['SAMPLE_TYPE']})"
        out.append(label)
        if xray.get("SOURCE"):
            out.append(xray["SOURCE"])
    return out


def variables_measured(path):
    """Read off which analyses actually exist beside this README.

    An experiment's variable is also stated by its kind's template, and the
    compaction in ``expsim_metadata.records`` drops it here when the two agree.
    """
    directory = Path(path).resolve().parent
    found = [label for name, label in VARIABLES.items() if (directory / name).is_file()]
    if any(directory.glob("*OrderParameters.json")):
        found.append("C-H bond order parameter")
    if any(directory.glob("*_FormFactor.json")):
        found.append("X-ray scattering form factor")

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
    """Return what this dataset is part of.

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


def enrich(block: dict, readme, path, kind, dois, names, inherited=None) -> dict:
    """Add every Bioschemas property derivable from local data.

    ``inherited`` is the block the record inherits from its templates. Nothing
    of it is written here; it only keeps the record's own keywords from
    repeating an inherited one.
    """
    subjects = block.pop("_subjects", [])
    # Composed, never fetched: a registry title names the paper or deposition, and
    # hundreds of records share one of those.
    block["name"] = compose_name(readme, path, kind, block)
    block["description"] = compose_description(readme, path, kind, dois, names)
    if dois.cited:
        # `identifier` is deliberately not written: the web frontend derives it
        # from the record's own DOI field. sameAs keeps the resolvable link, and
        # it names the cited DOI rather than the one that was looked up: what
        # this dataset is a copy of is the data, not the paper describing it.
        block["sameAs"] = f"https://doi.org/{dois.cited}"
    if kind == "simulations" and readme.get("SYSTEM"):
        block["alternateName"] = readme["SYSTEM"]

    block["keywords"] = dedupe_keywords(
        composition_keywords(readme, kind, names) + split_subjects(subjects),
        (inherited or {}).get("keywords") or (),
    )
    block["measurementTechnique"] = measurement_technique(readme, path, kind)
    block["variableMeasured"] = variables_measured(path)
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
