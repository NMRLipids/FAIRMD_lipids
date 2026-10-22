"""
Inheritance of ``bioschema_properties`` from the template of a record's kind.

Much of a record's Bioschemas block follows from the kind of record alone: every
simulation carries the same EDAM keywords, every order-parameter experiment the
same measured variable and the databank's own licence. Those values are written
once, in a ``bioschema_properties.yaml`` template in the kind's folder --
``Simulations``, ``experiments/OrderParameters``, ``experiments/FormFactors`` --
and inherited by every ``README.yaml`` of that kind. A README stores only what
is specific to its record.

The *resolved* block has two layers and no more: the kind's template, and the
record's own block applied on top of it. No folder between the two is read.
Applying the record on top of the template is a recursive dict update, with two
departures:

- lists are *extended*: inherited items first, then the record's own, without
  duplicates. A plain update would let a record's own ``keywords`` -- its
  molecule names -- replace the inherited ontology terms;
- the properties in :data:`REPLACED_WHOLE` are value objects or lists belonging
  to one record and are never merged into: the record's value replaces the
  inherited one whole. Merged key by key, a record's own ``license`` would keep
  the inherited licence's ``name`` and ``sameAs``.

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
"""Name of the template file a kind's folder may hold."""

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

    :param parent: the inherited block, a kind's template; its :data:`NOT_INHERITED` keys are dropped.
    :param child: the record's own block.
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


def template_path(template_dir) -> Path | None:
    """The template a kind's folder holds, or None when it holds none.

    :param template_dir: the kind's folder, e.g. ``Simulations`` or ``experiments/OrderParameters``.
    """
    path = Path(template_dir) / TEMPLATE_FILENAME
    return path if path.is_file() else None


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


def inherited_block(template_dir) -> dict:
    """The block every record of a kind inherits: its template, identity keys dropped.

    :param template_dir: the kind's folder, e.g. ``Simulations`` or ``experiments/OrderParameters``.
    :return: a new dict, empty when the folder holds no template.
    """
    path = template_path(template_dir)
    if path is None:
        return {}
    return {key: value for key, value in load_template(path).items() if key not in NOT_INHERITED}


def resolve_block(block: dict | None, template_dir) -> dict:
    """The record's complete Bioschemas block: its kind's template, updated by its own.

    :param block: the ``bioschema_properties`` block stored in the README.
    :param template_dir: the kind's folder, e.g. ``Simulations`` or ``experiments/OrderParameters``.
    :return: a new dict; the stored block is not modified.
    """
    return merge_blocks(inherited_block(template_dir), block)
