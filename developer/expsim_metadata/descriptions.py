"""
Composing the ``name`` and ``description`` of a simulation or experiment record.

Both are written from the record's own fields rather than fetched, because a
registry title names the paper or the Zenodo deposition and hundreds of records
share one of those. ``autocomplete_expsim_metadata`` calls in here for the two
composed properties and for nothing else.

A description is three sentences -- what the system is, how it was measured,
where it came from -- and there is one function per kind of analysis:

- :func:`simulation_description` -- molecular dynamics trajectories
- :func:`nmr_description`        -- order parameters measured by NMR
- :func:`xray_description`       -- form factors measured by scattering

The three share the sentences that do not depend on the technique, so the parts
that differ -- the instrument, the sample, whether a value was digitised out of
a paper -- are the only thing each function spells out. :func:`compose_description`
picks between them, which is the only entry point the caller needs.

Every value read out of the record comes from ``expsim_metadata.fields``; this module
only arranges those values into sentences.
"""

from pathlib import Path

from .constants import NULLISH
from .fields import chmo_for_method, composition_items, experiment_kind, solution_ids
from .helpers import clean_text, iso_date, number

# ---------------------------------------------------------------------------
# Phrases the name and the description are built from
# ---------------------------------------------------------------------------


def format_ratio(items) -> str:
    """Give human-readable composition: ``POPC`` for one component, ``POPC/POPE (95:5)`` for a mixture."""
    if not items:
        return "lipid"
    if len(items) == 1:
        return items[0][0]
    total = sum(amount for _, amount in items) or 1.0
    shares = ":".join(f"{100 * amount / total:.0f}" for _, amount in items)
    return "/".join(mol for mol, _ in items) + f" ({shares})"


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


def _bilayer(items, names, counted=False):
    """The component list every description opens with.

    Simulations state molecule counts, experiments molar fractions, so the count
    is written only where one was recorded. A single-component bilayer carries no
    share: ``100 mol%`` says nothing the list does not already say.
    """
    if not items:
        return "unrecorded composition"
    total = sum(amount for _, amount in items) or 1.0
    return _listed([
        _component(mol, names,
                   count=amount if counted else None,
                   share=None if len(items) == 1 else 100 * amount / total)
        for mol, amount in items
    ])


def _sentences(*parts):
    return " ".join(f"{part}." for part in parts)


# ---------------------------------------------------------------------------
# The composed name
# ---------------------------------------------------------------------------


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

    Simulations carry none: their databank ``ID`` is assigned only after merge,
    and the website's JSON-LD identifies each one by its page rather than by its
    name. Experiments use the first author and year -- which distinguishes the
    systems that two different groups measured independently -- falling back to
    the record's own directory when there is no publication to name.
    """
    if kind == "simulations":
        return None

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
    """Compose the record's title, pasted together from the record's own values."""
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


# ---------------------------------------------------------------------------
# The composed description, one function per kind of analysis
# ---------------------------------------------------------------------------


def simulation_description(readme, dois, names):
    """Describe a molecular dynamics trajectory.

    The middle sentence names what produced the trajectory -- force field,
    engine, length, size. A record that says none of that gets ``Trajectory
    run``: claiming a method nothing in the record supports would be worse than
    saying little.
    """
    system = ("Molecular dynamics trajectory of a lipid bilayer of "
              + _bilayer(composition_items(readme, "simulations"), names, counted=True))
    temperature = number(readme.get("TEMPERATURE"))
    if temperature:
        system += f" at {temperature} K"

    force_field = clean_text(readme.get("FF"))
    engine = clean_text(readme.get("SOFTWARE"))
    if engine:
        version = clean_text(readme.get("SOFTWARE_VERSION"))
        engine = engine.upper() + (f" {version}" if version else "")
    method = "Simulated"
    if force_field:
        method += f" with {force_field}"
    if engine:
        method += f" in {engine}"
    if not force_field and not engine:
        # No record says what produced it, so do not claim a method.
        method = "Trajectory run"
    if readme.get("TRJLENGTH"):
        method += f" for {float(readme['TRJLENGTH']) / 1000:.0f} ns"
    if readme.get("NUMBER_OF_ATOMS"):
        method += f" ({readme['NUMBER_OF_ATOMS']} atoms)"

    origin = "Part of the NMRlipids Databank"
    if dois.cited:
        origin += f" and available from https://doi.org/{dois.cited}"
    return _sentences(system, method, origin)


def nmr_description(readme, path, dois, names):
    """Describe order parameters measured by NMR.

    The method sentence names the CHMO term the record maps onto, then the two
    things that tell two NMR measurements of one system apart: the pulse
    sequence and the spectrometer.
    """
    method = f"Measured by {chmo_for_method(readme, path)[1]}"
    nmr = readme.get("NMR") or {}
    if nmr.get("METHOD"):
        method += f" ({clean_text(nmr['METHOD'])})"
    if nmr.get("INSTRUMENT"):
        method += f" on a {clean_text(nmr['INSTRUMENT'])}"
    return _sentences(_sample_sentence(readme, path, names), method, _origin_sentence(dois))


def xray_description(readme, path, dois, names):
    """Describe a form factor measured by small-angle X-ray scattering.

    Scattering records state the sample geometry and the beamline instead of a
    pulse sequence, and those are what separate two measurements of one system.
    """
    method = f"Measured by {chmo_for_method(readme, path)[1]}"
    xray = readme.get("XRAY") or {}
    if xray.get("SAMPLE_TYPE"):
        method += f" on {clean_text(xray['SAMPLE_TYPE'])} samples"
    if xray.get("SOURCE"):
        method += f" at {clean_text(xray['SOURCE'])}"
    return _sentences(_sample_sentence(readme, path, names), method, _origin_sentence(dois))


def _sample_sentence(readme, path, names):
    """What an experiment was measured on: composition, temperature, conditions.

    Shared by every experimental technique, since the sample is a property of the
    record and not of the instrument pointed at it.
    """
    sentence = (f"Experimental {measured_quantity(path)} for a lipid bilayer of "
                + _bilayer(composition_items(readme, "experiments"), names))
    temperature = number(readme.get("TEMPERATURE"))
    if temperature:
        sentence += f" at {temperature} K"
    hydration = hydration_phrase(readme)
    if hydration == "full hydration":
        sentence += ", fully hydrated"
    elif hydration:
        sentence += f", hydrated to {hydration}"

    conditions = []
    ions = solution_ids(readme)
    if ions:
        conditions.append("ions " + _listed([_named(ion, names) for ion in ions]))
    additives = sorted((readme.get("ADDITIONAL_MOLECULES") or {}).keys())
    if additives:
        conditions.append("additives " + _listed(additives))
    ph = readme.get("PH")
    if ph is not None and str(ph).strip().lower() not in NULLISH | {"unknown"}:
        how = clean_text(readme.get("PH_METHOD"))
        conditions.append(f"pH {number(ph)}"
                          + (f" ({how})" if how and how.lower() != "unknown" else ""))
    if conditions:
        sentence += ", with " + "; ".join(conditions)
    return sentence


def _origin_sentence(dois):
    """Where an experimental record's values came from.

    The cited DOI, so this sentence, ``sameAs`` and ``citation`` all name one
    source rather than three.
    """
    if dois.cited:
        return f"Values digitised into the NMRlipids Databank from https://doi.org/{dois.cited}"
    return "Unpublished data contributed to the NMRlipids Databank"


def compose_description(readme, path, kind, dois, names):
    """Three sentences: what the system is, how it was measured, where it came from.

    The record's own directory decides which analysis it holds, exactly as
    ``chmo_for_method`` decides the ontology term: an OrderParameters record with
    no ``NMR:`` block is still an NMR measurement, and around twenty of them
    carry none.
    """
    if kind == "simulations":
        return simulation_description(readme, dois, names)
    if experiment_kind(path) == "xray":
        return xray_description(readme, path, dois, names)
    return nmr_description(readme, path, dois, names)
