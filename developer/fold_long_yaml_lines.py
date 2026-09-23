"""
Fold ``README.yaml`` values that run past ``FOLD_WIDTH`` columns into ``>`` blocks.

This is the hand-written counterpart of ``BlockDumper`` in
``expsim_metadata.records``: the generated ``bioschema_properties`` block is
folded when the autocomplete writes it, and this script folds the rest of a
README in the same way. BilayerData's ``LintReadmeLineLength`` workflow runs it
on the READMEs a pull request changes and posts the result as review
suggestions; run over the whole databank it is also the one-off cleanup of the
lines written before that check existed.

Like the autocomplete, it never re-serialises a file. Only the lines holding an
overlong value are rewritten, and an edit is kept only when the file still parses
to exactly the same data -- so comments, quoting and layout elsewhere survive,
and a value that cannot be folded without changing it is left alone:

- a single-line value -- plain or quoted -- becomes a ``>-`` block, with a
  trailing comment moved onto the indicator line;
- an existing ``>`` block is re-wrapped, keeping its indicator;
- literal ``|`` blocks, values written over several lines, and values with no
  space in them are reported and left as they are.

Usage::

    python developer/fold_long_yaml_lines.py [--check] README.yaml [...]

``--check`` writes nothing and exits 1 when a file would change.
"""

import argparse
import importlib.util
import sys
from pathlib import Path

import yaml


def _load_encoder():
    """Load the formatter module without importing ``fairmd.lipids``."""
    module_path = (
        Path(__file__).resolve().parents[1]
        / "src"
        / "fairmd"
        / "lipids"
        / "auxiliary"
        / "yaml_format.py"
    )
    spec = importlib.util.spec_from_file_location("fairmd_yaml_format", module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load YAML formatter from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.encode_canonical_yaml


encode_canonical_yaml = _load_encoder()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+", type=Path, help="YAML files to fold in place")
    parser.add_argument("--check", action="store_true", help="write nothing; exit 1 if a file would change")
    args = parser.parse_args()

    changed = 0
    for path in args.files:
        original = path.read_text(encoding="utf-8")
        folded = encode_canonical_yaml(yaml.safe_load(original))
        if folded != original:
            changed += 1
            if not args.check:
                path.write_text(folded, encoding="utf-8")
            print(f"{'Would fold' if args.check else 'Folded'} long lines in {path}")
    return 1 if args.check and changed else 0


if __name__ == "__main__":
    sys.exit(main())
