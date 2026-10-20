"""
Inheritance of ``bioschema_properties`` from the folders above a record.

Much of a record's Bioschemas block follows from the kind of record alone: every
simulation carries the same EDAM keywords, every order-parameter experiment the
same measured variable, every experiment the databank's own licence. Those values
are written once, in a ``bioschema_properties.yaml`` template in a folder, and
inherited by every ``README.yaml`` below it. A README stores only what is
specific to its record.

The *resolved* block is the record's own block applied on top of the templates
of every folder from the kind base (``Simulations`` or ``experiments``) down to
the record's parent folder, outermost first. Applying one level on top of another
is a recursive dict update, with two departures:

- lists are *extended*: inherited items first, then the closer level's own,
  without duplicates. A plain update would let a record's own ``keywords`` -- its
  molecule names -- replace the inherited ontology terms;
- the properties in :data:`REPLACED_WHOLE` are value objects or lists belonging
  to one record and are never merged into: a closer value replaces the inherited
  one whole. Merged key by key, a record's own ``license`` would keep the
  inherited licence's ``name`` and ``sameAs``.

The properties in :data:`NOT_INHERITED` identify one record or record where its
metadata came from, and a template value for them is ignored.

Only the resolved block is a complete Bioschemas ``Dataset``. It is computed on
demand and never written back into a README.

.. note::
   ``developer/expsim_metadata`` loads this file by path, so the metadata
   generator runs without the package installed. It must import nothing from
   ``fairmd``: the standard library and PyYAML only.
"""

import copy
from functools import lru_cache
from pathlib import Path

import yaml

TEMPLATE_FILENAME = "bioschema_properties.yaml"
"""Name of the template file a folder may hold."""

NOT_INHERITED = frozenset({"_source", "name", "alternateName", "description", "sameAs"})
"""Top-level properties a template may not set: the identity and provenance of one record."""

REPLACED_WHOLE = frozenset(
    {"license", "articleLicense", "isPartOf", "creator", "citation", "distribution", "isBasedOn"},
)
"""Top-level properties inherited only as a whole: a closer value replaces, never merges."""

_DROP = object()


def _update(base, override):
    """Recursive update of ``base`` by ``override``, lists extended; inputs are not modified."""
    if isinstance(base, dict) and isinstance(override, dict):
        out = {}
        for key, value in override.items():
            out[key] = _update(base[key], value) if key in base else copy.deepcopy(value)
        for key, value in base.items():
            if key not in out:
                out[key] = copy.deepcopy(value)
        return out
    if isinstance(base, list) and isinstance(override, list):
        out = copy.deepcopy(base)
        for item in override:
            if item not in out:
                out.append(copy.deepcopy(item))
        return out
    return copy.deepcopy(override)


def merge_blocks(parent: dict | None, child: dict | None) -> dict:
    """Apply ``child`` on top of the inherited ``parent`` block.

    :param parent: the inherited block; its :data:`NOT_INHERITED` keys are dropped.
    :param child: the closer level, a template or the record's own block.
    :return: a new dict; neither input is modified.
    """
    parent = {key: value for key, value in (parent or {}).items() if key not in NOT_INHERITED}
    out = {}
    for key, value in (child or {}).items():
        if key in parent and key not in REPLACED_WHOLE:
            out[key] = _update(parent[key], value)
        else:
            out[key] = copy.deepcopy(value)
    for key, value in parent.items():
        if key not in out:
            out[key] = copy.deepcopy(value)
    return out


def _subtract(value, base):
    """What of ``value`` is not already stated by ``base``, or ``_DROP`` for nothing."""
    if isinstance(value, dict) and isinstance(base, dict):
        out = {}
        for key, item in value.items():
            rest = _subtract(item, base[key]) if key in base else item
            if rest is not _DROP:
                out[key] = rest
        return out or _DROP
    if isinstance(value, list) and isinstance(base, list):
        out = [item for item in value if item not in base]
        return out or _DROP
    return _DROP if value == base else value


def compact_block(block: dict, inherited: dict | None) -> dict:
    """Strip from ``block`` what it would inherit anyway; the inverse of :func:`merge_blocks`.

    ``merge_blocks(inherited, compact_block(block, inherited))`` equals
    ``merge_blocks(inherited, block)``: compacting loses nothing the resolved
    block states.

    :param block: a full record block.
    :param inherited: the block the record inherits, see :func:`inherited_block`.
    :return: a new dict holding only what the record has to state itself.
    """
    inherited = {key: value for key, value in (inherited or {}).items() if key not in NOT_INHERITED}
    out = {}
    for key, value in block.items():
        if key not in inherited:
            out[key] = copy.deepcopy(value)
        elif key in REPLACED_WHOLE:
            if value != inherited[key]:
                out[key] = copy.deepcopy(value)
        else:
            rest = _subtract(value, inherited[key])
            if rest is not _DROP:
                out[key] = copy.deepcopy(rest)
    return out


def template_paths(record_dir, base_dir) -> list[Path]:
    """The templates a record inherits from, outermost first.

    :param record_dir: the folder holding the record's ``README.yaml``.
    :param base_dir: the kind base, e.g. the ``Simulations`` or ``experiments`` folder.
    :return: existing template files from ``base_dir`` down to the record's parent
             folder; empty when the record is not below ``base_dir``.
    """
    record_dir = Path(record_dir).resolve()
    base_dir = Path(base_dir).resolve()
    try:
        relative = record_dir.relative_to(base_dir)
    except ValueError:
        return []
    if not relative.parts:
        return []
    folders = [base_dir]
    for part in relative.parts[:-1]:
        folders.append(folders[-1] / part)
    return [folder / TEMPLATE_FILENAME for folder in folders if (folder / TEMPLATE_FILENAME).is_file()]


@lru_cache(maxsize=None)
def _read_template(path: str, mtime_ns: int) -> dict:
    # mtime_ns is unused here: it is part of the cache key, so an edited template is re-read.
    with open(path, encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if data is None:
        return {}
    if not isinstance(data, dict):
        msg = f"{path}: a {TEMPLATE_FILENAME} template must be a mapping of Bioschemas properties"
        raise ValueError(msg)
    return data


def load_template(path) -> dict:
    """Read one template file, parsed once per modification and returned as a copy."""
    path = Path(path)
    return copy.deepcopy(_read_template(str(path), path.stat().st_mtime_ns))


def inherited_block(record_dir, base_dir) -> dict:
    """The block a record inherits: its templates merged outermost first.

    :param record_dir: the folder holding the record's ``README.yaml``.
    :param base_dir: the kind base, e.g. the ``Simulations`` or ``experiments`` folder.
    :return: a new dict, empty when no template applies.
    """
    inherited = {}
    for path in template_paths(record_dir, base_dir):
        template = {key: value for key, value in load_template(path).items() if key not in NOT_INHERITED}
        inherited = merge_blocks(inherited, template)
    return inherited


def resolve_block(block: dict | None, record_dir, base_dir) -> dict:
    """The record's complete Bioschemas block: what it inherits, updated by its own.

    :param block: the ``bioschema_properties`` block stored in the README.
    :param record_dir: the folder holding the record's ``README.yaml``.
    :param base_dir: the kind base, e.g. the ``Simulations`` or ``experiments`` folder.
    :return: a new dict; the stored block is not modified.
    """
    return merge_blocks(inherited_block(record_dir, base_dir), block)
