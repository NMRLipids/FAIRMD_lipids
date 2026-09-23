"""
A library for formatting YAML files according to the project's style guide.

Exports one main function, :func:`encode_canonical_yaml`, which serializes a Python object to YAML
and folds any values that run past ``FOLD_WIDTH`` columns into ``>`` blocks.
"""

import re
import textwrap

import yaml

FOLD_WIDTH = 100
FOLDABLE = re.compile(r"\S+( \S+)+")

STR_TAG = "tag:yaml.org,2002:str"


def _value_scalars(node, parent_column=None):
    """Every block-context scalar value with the column its folded lines go to.

    That column is two past the key's, or two past a sequence item's dash, which
    is where ``BlockDumper`` puts them too. Keys and flow collections are skipped.
    """
    if isinstance(node, yaml.MappingNode) and not node.flow_style:
        for key, value in node.value:
            yield from _value_scalars(value, key.start_mark.column)
    elif isinstance(node, yaml.SequenceNode) and not node.flow_style:
        for item in node.value:
            yield from _value_scalars(item, item.start_mark.column - 2)
    elif isinstance(node, yaml.ScalarNode) and parent_column is not None and node.tag == STR_TAG:
        yield node, parent_column + 2


def _wrap(text, indent):
    return [" " * indent + line + "\n"
            for line in textwrap.wrap(text, FOLD_WIDTH - indent,
                                      break_long_words=False, break_on_hyphens=False)]


def _too_long(line):
    return len(line.rstrip("\r\n")) > FOLD_WIDTH


def _fold_single_line(node, indent, lines):
    """``key: long value`` -> ``key: >-`` followed by the folded value."""
    start, end = node.start_mark, node.end_mark
    line = lines[start.line]
    if not _too_long(line):
        return None
    if start.line != end.line:
        return "written over several lines"
    if not FOLDABLE.fullmatch(node.value):
        return "no space to fold at"
    comment = line[end.column:].strip()
    if comment and not comment.startswith("#"):
        return "unexpected text after the value"
    head = line[:start.column] + ">-" + (f"  {comment}" if comment else "") + "\n"
    return start.line, start.line + 1, [head] + _wrap(node.value, indent)


def _rewrap_folded(node, lines):
    """Re-wrap the content of an existing ``>`` block, keeping its indicator line."""
    first, stop = node.start_mark.line + 1, node.end_mark.line
    while stop > first and not lines[stop - 1].strip():
        stop -= 1  # blank lines after the text stay where they are
    content = lines[first:stop]
    if not any(_too_long(line) for line in content):
        return None
    paragraphs = node.value.rstrip("\n").split("\n")
    if not all(FOLDABLE.fullmatch(p) for p in paragraphs):
        return "cannot be re-wrapped without changing it"
    indent = len(content[0]) - len(content[0].lstrip(" "))
    new = []
    for paragraph in paragraphs:
        if new:
            new.append("\n")  # one blank line is what separates two paragraphs
        new += _wrap(paragraph, indent)
    return first, stop, new


def _fold_file(text):
    """The text with overlong values folded, and ``(line, reason)`` for those left."""
    root = yaml.compose(text, Loader=yaml.SafeLoader)
    if root is None:
        return text, []
    data = yaml.safe_load(text)
    lines = text.splitlines(keepends=True)

    edits, left = [], []
    for node, indent in _value_scalars(root):
        if node.style == ">":
            result = _rewrap_folded(node, lines)
        elif node.style == "|":
            result = "literal block" if any(_too_long(line) for line in lines[
                node.start_mark.line + 1:node.end_mark.line]) else None
        else:
            result = _fold_single_line(node, indent, lines)
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


def encode_canonical_yaml(data) -> str:
    """Serialize data using the project YAML formatting rules."""
    text = yaml.dump(data, sort_keys=False, allow_unicode=True,
                     default_flow_style=False, width=4096)
    return _fold_file(text)[0]
