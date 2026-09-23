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
import sys
import textwrap
from pathlib import Path

import yaml

try:  # run as a script: developer/ is on sys.path
    from expsim_metadata.constants import FOLD_WIDTH, FOLDABLE
except ImportError:  # imported as developer.fold_long_yaml_lines, e.g. by autodoc
    from developer.expsim_metadata.constants import FOLD_WIDTH, FOLDABLE

STR_TAG = "tag:yaml.org,2002:str"


def value_scalars(node, parent_column=None):
    """Every block-context scalar value with the column its folded lines go to.

    That column is two past the key's, or two past a sequence item's dash, which
    is where ``BlockDumper`` puts them too. Keys and flow collections are skipped.
    """
    if isinstance(node, yaml.MappingNode) and not node.flow_style:
        for key, value in node.value:
            yield from value_scalars(value, key.start_mark.column)
    elif isinstance(node, yaml.SequenceNode) and not node.flow_style:
        for item in node.value:
            yield from value_scalars(item, item.start_mark.column - 2)
    elif isinstance(node, yaml.ScalarNode) and parent_column is not None and node.tag == STR_TAG:
        yield node, parent_column + 2


def wrap(text, indent):
    return [" " * indent + line + "\n"
            for line in textwrap.wrap(text, FOLD_WIDTH - indent,
                                      break_long_words=False, break_on_hyphens=False)]


def too_long(line):
    return len(line.rstrip("\r\n")) > FOLD_WIDTH


def fold_single_line(node, indent, lines):
    """``key: long value`` -> ``key: >-`` followed by the folded value."""
    start, end = node.start_mark, node.end_mark
    line = lines[start.line]
    if not too_long(line):
        return None
    if start.line != end.line:
        return "written over several lines"
    if not FOLDABLE.fullmatch(node.value):
        return "no space to fold at"
    comment = line[end.column:].strip()
    if comment and not comment.startswith("#"):
        return "unexpected text after the value"
    head = line[:start.column] + ">-" + (f"  {comment}" if comment else "") + "\n"
    return start.line, start.line + 1, [head] + wrap(node.value, indent)


def rewrap_folded(node, lines):
    """Re-wrap the content of an existing ``>`` block, keeping its indicator line."""
    first, stop = node.start_mark.line + 1, node.end_mark.line
    while stop > first and not lines[stop - 1].strip():
        stop -= 1  # blank lines after the text stay where they are
    content = lines[first:stop]
    if not any(too_long(line) for line in content):
        return None
    paragraphs = node.value.rstrip("\n").split("\n")
    if not all(FOLDABLE.fullmatch(p) for p in paragraphs):
        return "cannot be re-wrapped without changing it"
    indent = len(content[0]) - len(content[0].lstrip(" "))
    new = []
    for paragraph in paragraphs:
        if new:
            new.append("\n")  # one blank line is what separates two paragraphs
        new += wrap(paragraph, indent)
    return first, stop, new


def fold_file(text):
    """The text with overlong values folded, and ``(line, reason)`` for those left."""
    root = yaml.compose(text, Loader=yaml.SafeLoader)
    if root is None:
        return text, []
    data = yaml.safe_load(text)
    lines = text.splitlines(keepends=True)

    edits, left = [], []
    for node, indent in value_scalars(root):
        if node.style == ">":
            result = rewrap_folded(node, lines)
        elif node.style == "|":
            result = "literal block" if any(too_long(line) for line in lines[
                node.start_mark.line + 1:node.end_mark.line]) else None
        else:
            result = fold_single_line(node, indent, lines)
        if isinstance(result, str):
            left.append((node.start_mark.line + 1, result))
        elif result:
            edits.append(result)

    # Bottom up, so an edit never moves the lines a later one refers to. Each is
    # kept only if the file still reads back as the same data.
    for first, stop, new in sorted(edits, reverse=True):
        candidate = lines[:first] + new + lines[stop:]
        if yaml.safe_load("".join(candidate)) == data:
            lines = candidate
        else:
            left.append((first + 1, "folding would change the value"))
    return "".join(lines), sorted(left)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+", type=Path, help="YAML files to fold in place")
    parser.add_argument("--check", action="store_true", help="write nothing; exit 1 if a file would change")
    args = parser.parse_args()

    changed = 0
    for path in args.files:
        original = path.read_text(encoding="utf-8")
        folded, left = fold_file(original)
        for line, reason in left:
            print(f"{path}:{line}: left as is ({reason})")
        if folded != original:
            changed += 1
            if not args.check:
                path.write_text(folded, encoding="utf-8")
            print(f"{'Would fold' if args.check else 'Folded'} long lines in {path}")
    return 1 if args.check and changed else 0


if __name__ == "__main__":
    sys.exit(main())
