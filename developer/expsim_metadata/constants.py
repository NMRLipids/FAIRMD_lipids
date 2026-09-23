"""
Constants shared by the simulation and experiment autocomplete modules.

Every value here is a fact about the outside world -- an endpoint, an ontology
term, a licence, the shape of the block written into a ``README.yaml`` -- rather
than about one record. They live together so that a term can be checked or
updated in one place, and so that the modules using them stay free of literals.

Nothing in this module is computed and nothing imports from its siblings.
"""

import os
import re
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

# The terms and the databank's own names are data, not logic, and live in a file
# a curator can read: see constants.yaml for what was chosen for what, and why.
CONSTANTS_FILE = Path(__file__).with_name("constants.yaml")
_DATA = yaml.safe_load(CONSTANTS_FILE.read_text(encoding="utf-8"))


def _term(*keys):
    """One ``(termCode, name)`` pair from constants.yaml."""
    entry = _DATA
    for key in keys:
        entry = entry[key]
    return tuple(entry)


def _unit(entry):
    """One ``(uri, symbol)`` pair from the ``units`` section of constants.yaml.

    A unit's term is written ``<set>/<code>`` against the sets in the same file.
    A quantity with no term keeps its plain-text symbol alone.
    """
    term = entry.get("term")
    if not term:
        return (None, entry["symbol"])
    term_set, code = term.split("/", 1)
    return (f"{_DATA['sets'][term_set]}/{code}", entry["symbol"])


def _terms(*keys):
    """A list of ``(termCode, name)`` pairs from constants.yaml."""
    entry = _DATA
    for key in keys:
        entry = entry[key]
    return [tuple(term) for term in entry]


# The two sets the generated terms name themselves against; the QUDT set is
# only ever reached through a unit, so it stays in constants.yaml.
EDAM_SET = _DATA["sets"]["edam"]
OBO_SET = _DATA["sets"]["obo"]

EDAM_SIMULATION = _terms("edam_keywords", "simulation")
EDAM_NMR = _terms("edam_keywords", "nmr")
EDAM_XRAY = _terms("edam_keywords", "xray")

EDAM_MD_OP = _term("edam_techniques", "md_operation")
EDAM_NMR_TOPIC = _term("edam_techniques", "nmr_topic")
EDAM_XRAY_TOPIC = _term("edam_techniques", "xray_topic")

CHMO_SAXS = _term("chmo", "saxs")
CHMO_NMR = _term("chmo", "nmr")
CHMO_SSNMR = _term("chmo", "ssnmr")
CHMO_PDLF = _term("chmo", "pdlf")
CHMO_PDLF_R = _term("chmo", "pdlf_r")

# The units of the measured quantities, as ``(uri, symbol)``; which term was
# picked for what, and why one quantity has none, is in constants.yaml.
UNITS = {name: _unit(entry) for name, entry in _DATA["units"].items()}

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

# Which analysis outputs correspond to which measured quantity, and the media
# types of the trajectory formats; both are the databank's own names, listed in
# the second half of constants.yaml.
VARIABLES = dict(_DATA["variables"])
TRAJECTORY_FORMATS = dict(_DATA["trajectory_formats"])

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

# Prose in the generated block that would run past this column -- composed names
# and descriptions, registry titles -- is folded into a ``>-`` block scalar.
FOLD_WIDTH = 100
# Words separated by single spaces: the only text a fold reads back unchanged,
# since the parser turns each line break back into exactly one space.
FOLDABLE = re.compile(r"\S+( \S+)+")

DOI_RE = re.compile(r"10\.\d{4,9}/[^\s;,\"']+")
ORCID_RE = re.compile(r"(\d{4}-\d{4}-\d{4}-\d{3}[\dX])", re.IGNORECASE)
DATE_RE = re.compile(r"^\d{4}(-\d{2}(-\d{2})?)?$")
TRAILING = ".,;:)]}’\"'"


# Experiment README keys that are no longer accepted, mapped to the spelling
# that replaces them. They are still read, so a record that has not been
# migrated keeps its citation, but every read is reported: experiment_schema.json
# declares none of them while setting additionalProperties: false, so a record
# still carrying one does not validate. Renaming the key is an edit to
# hand-written content, which this tool does not make -- it is done in the data
# repository.
#
# DOI is deprecated for experiments only: simulation_schema.json declares DOI
# for a simulation, where it names the Zenodo deposition and is the right key.
DEPRECATED_EXPERIMENT_KEYS = {
    "DOI": "ARTICLE_DOI, or DATA_DOI where the value is a data deposition",
}


# Water is in every system and says nothing about it; ions do carry meaning.
KEYWORD_SKIP = frozenset({"SOL"})

# Values that mean "this record does not say", written half a dozen ways.
NULLISH = frozenset({"null", "none", "na", "n/a", "??", "?", ""})

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
# A simulation's DATE is not declared in simulation_schema.json at all, so it
# is left as YAML parses it rather than rewritten for tidiness.
DATE_FIELDS = {
    "simulations": {"DATEOFRUNNING": True, "FF_DATE": False},
    "experiments": {"DATE": False},
}
