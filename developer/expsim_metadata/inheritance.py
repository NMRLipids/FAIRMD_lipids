"""
The inheritance rules, shared with the package rather than copied.

What a record inherits from the ``bioschema_properties.yaml`` templates above it,
and how a block is compacted against them, is decided in one place:
``src/fairmd/lipids/bioschema_inheritance.py``. The package resolves blocks with
it when a record is loaded, and this tool compacts them with it when one is
written, so the two can never disagree.

The file is loaded by path rather than imported as ``fairmd.lipids...``: the
tool runs from a plain checkout, without the package installed. That file
imports nothing from ``fairmd`` for exactly this reason.
"""

import importlib.util
import sys
from pathlib import Path

from .fields import kind_base_of

SOURCE = Path(__file__).resolve().parents[2] / "src" / "fairmd" / "lipids" / "bioschema_inheritance.py"
_MODULE_NAME = "_expsim_bioschema_inheritance"


def _load():
    if _MODULE_NAME in sys.modules:
        return sys.modules[_MODULE_NAME]
    spec = importlib.util.spec_from_file_location(_MODULE_NAME, SOURCE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[_MODULE_NAME] = module
    spec.loader.exec_module(module)
    return module


_rules = _load()

TEMPLATE_FILENAME = _rules.TEMPLATE_FILENAME
NOT_INHERITED = _rules.NOT_INHERITED
REPLACED_WHOLE = _rules.REPLACED_WHOLE
compact_block = _rules.compact_block
load_template = _rules.load_template
merge_blocks = _rules.merge_blocks


def is_template(path):
    return Path(path).name == TEMPLATE_FILENAME


def templates_for(path):
    """The templates a ``README.yaml`` inherits from, outermost first."""
    record_dir = Path(path).resolve().parent
    return _rules.template_paths(record_dir, kind_base_of(path))


def inherited_for(path):
    """The block a ``README.yaml`` inherits from the templates above it."""
    record_dir = Path(path).resolve().parent
    return _rules.inherited_block(record_dir, kind_base_of(path))


def resolved_for(path, block):
    """A ``README.yaml``'s complete block: what it inherits, updated by its own."""
    return merge_blocks(inherited_for(path), block)
