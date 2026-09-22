"""
Reading facts out of a simulation or experiment ``README.yaml``.

The values here are what a record states about itself -- which kind of analysis
it holds, which DOIs it carries, what the bilayer is made of, which method was
pointed at it -- normalised so that two records saying the same thing in
different words give the same value.

Nothing in this module composes prose or Bioschemas properties: it is the layer
``expsim_metadata.descriptions`` and ``expsim_metadata.bioschema`` both read the record through.
"""

from collections import namedtuple
from pathlib import Path

import yaml

from .constants import (
    CHMO_NMR,
    CHMO_PDLF,
    CHMO_PDLF_R,
    CHMO_SAXS,
    CHMO_SSNMR,
    DEPRECATED_EXPERIMENT_KEYS,
    KEYWORD_SKIP,
)
from .helpers import normalize_doi

# ---------------------------------------------------------------------------
# Where a record sits and what it points at
# ---------------------------------------------------------------------------


SIMULATIONS_DIR = "Simulations"


def _is_simulations_dir(name):
    """``Simulations``, and the ``Simulations.1`` variants the toy data uses.

    The databank keeps its trajectories in one ``Simulations`` folder, but the
    test data shipped with the package splits them across ``Simulations.1``,
    ``Simulations.2`` and ``Simulations.AddData``. Matching only the bare name
    files every one of those records as an experiment, which then picks the
    wrong DOI roles, the wrong description and the wrong schema.
    """
    return name == SIMULATIONS_DIR or name.startswith(SIMULATIONS_DIR + ".")


def record_kind(path):
    parts = Path(path).resolve().parts
    return "simulations" if any(_is_simulations_dir(part) for part in parts) else "experiments"


def experiment_kind(path):
    return "xray" if "FormFactors" in Path(path).resolve().parts else "nmr"


def data_root_of(path):
    """Walk up from a README.yaml to the databank root."""
    for parent in Path(path).resolve().parents:
        if not (parent / "Molecules").is_dir():
            continue
        if any(child.is_dir() and _is_simulations_dir(child.name)
               for child in parent.glob(f"{SIMULATIONS_DIR}*")):
            return parent
    # Both layouts put README.yaml five levels below the root.
    return Path(path).resolve().parents[5]


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


DOI_ROLES = ("lookup", "cited", "article", "deposition")
Dois = namedtuple("Dois", DOI_ROLES)


def record_dois(readme, kind):
    """Get record's DOIs, split by the role each one plays.

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

    ``DOI`` is deprecated for an experiment, in favour of ``ARTICLE_DOI`` and
    ``DATA_DOI``. It is still read so that the records carrying it keep their
    citation, but it is no longer a spelling a record may use: see
    ``deprecated_keys`` below, which is what reports it.
    """
    if kind == "simulations":
        # A simulation has one DOI, the Zenodo deposition holding the
        # trajectory. Nothing is cited from it: the deposition *is* this record,
        # so it belongs in isPartOf and distribution, not in a citation list.
        deposition = normalize_doi(readme.get("DOI"))
        return Dois(lookup=deposition, cited=deposition,
                    article=None, deposition=deposition)

    # readme.get("DOI") is the deprecated spelling, kept only until the records
    # carrying it are migrated.
    article = normalize_doi(readme.get("ARTICLE_DOI") or readme.get("DOI"))
    deposition = normalize_doi(readme.get("DATA_DOI"))
    return Dois(lookup=article or deposition, cited=deposition or article,
                article=article, deposition=deposition)


def deprecated_keys(readme, kind):
    """Deprecated spellings this record still uses, as ``{key: replacement}``.

    Reported rather than rewritten. ``experiment_schema.json`` declares none of
    these while setting ``additionalProperties: false``, so a record still
    carrying one does not validate, and the fix is to rename the key in the
    record -- an edit to hand-written content, which this tool does not make.
    What it does instead is name the key every time it reads one, so a run over
    the databank lists exactly which records are still to be migrated.

    A simulation has none: ``DOI`` is the deprecated *experiment* spelling, and
    the same key on a simulation is the Zenodo deposition that
    ``readme_yaml_schema.json`` declares.
    """
    if kind != "experiments":
        return {}
    return {key: replacement
            for key, replacement in DEPRECATED_EXPERIMENT_KEYS.items()
            if readme.get(key) is not None}


# ---------------------------------------------------------------------------
# What the record says it contains
# ---------------------------------------------------------------------------


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


def composition_items(readme, kind: str) -> list:
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


def solution_ids(readme) -> list:
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


def chmo_for_method(readme, path: str) -> str:
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
