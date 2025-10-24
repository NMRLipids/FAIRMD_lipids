#!/usr/bin/env python3

#
# Compact JSON from @jannismain is used here.
# https://gist.github.com/jannismain/e96666ca4f059c3e5bc28abb711b5c92
#

import json
import math

import numpy as np


class CompactJSONEncoder(json.JSONEncoder):
    """A JSON Encoder that puts small containers on single lines."""

    CONTAINER_TYPES = (list, tuple, dict)
    """Container datatypes include primitives or other containers."""

    MAX_WIDTH = 150
    """Maximum width of a container that might be put on a single line."""

    MAX_ITEMS = 10
    """Maximum number of items in container that might be put on single line."""

    def __init__(self, *args, **kwargs):
        # using this class without indentation is pointless
        if kwargs.get("indent") is None:
            kwargs["indent"] = 4
        super().__init__(*args, **kwargs)
        self.indentation_level = 0

    def encode(self, o):
        """Encode JSON object *o* with respect to single line lists."""
        if isinstance(o, (list, tuple)):
            return self._encode_list(o)
        if isinstance(o, dict):
            return self._encode_object(o)
        if isinstance(o, (float, np.floating)):  # Use scientific notation for floats
            if math.isnan(o):
                # !very important. Otherwise it will put "nan" with lowercase,
                # which is not recognized by standard json DEcoder
                return "NaN"
            return format(o, "g")
        return json.dumps(
            o,
            skipkeys=self.skipkeys,
            ensure_ascii=self.ensure_ascii,
            check_circular=self.check_circular,
            allow_nan=self.allow_nan,
            sort_keys=self.sort_keys,
            indent=self.indent,
            separators=(self.item_separator, self.key_separator),
            default=self.default if hasattr(self, "default") else None,
        )

    def _encode_list(self, o):
        if self._put_on_single_line(o):
            return "[" + ", ".join(self.encode(el) for el in o) + "]"
        self.indentation_level += 1
        output = [self.indent_str + self.encode(el) for el in o]
        self.indentation_level -= 1
        return "[\n" + ",\n".join(output) + "\n" + self.indent_str + "]"

    def _encode_object(self, o):
        if not o:
            return "{}"

        # ensure keys are converted to strings
        o = {str(k) if k is not None else "null": v for k, v in o.items()}

        if self.sort_keys:
            o = dict(sorted(o.items(), key=lambda x: x[0]))

        if self._put_on_single_line(o):
            return "{ " + ", ".join(f"{json.dumps(k)}: {self.encode(el)}" for k, el in o.items()) + " }"

        self.indentation_level += 1
        output = [f"{self.indent_str}{json.dumps(k)}: {self.encode(v)}" for k, v in o.items()]
        self.indentation_level -= 1

        return "{\n" + ",\n".join(output) + "\n" + self.indent_str + "}"

    def iterencode(self, o, **kwargs):
        """Required to also work with `json.dump`."""
        return self.encode(o)

    def _put_on_single_line(self, o):
        return self._primitives_only(o) and len(o) <= self.MAX_ITEMS and len(str(o)) - 2 <= self.MAX_WIDTH

    def _primitives_only(self, o: list | tuple | dict):
        if isinstance(o, (list, tuple)):
            return not any(isinstance(el, self.CONTAINER_TYPES) for el in o)
        if isinstance(o, dict):
            return not any(isinstance(el, self.CONTAINER_TYPES) for el in o.values())

    @property
    def indent_str(self) -> str:
        if isinstance(self.indent, int):
            return " " * (self.indentation_level * self.indent)
        if isinstance(self.indent, str):
            return self.indentation_level * self.indent
        raise ValueError(
            f"indent must either be of type int or str (is: {type(self.indent)})",
        )


if __name__ == "__main__":
    import sys

    if "--example" in sys.argv:
        data = {
            "compact_object": {"first": "element", "second": 2},
            "compact_list": ["first", "second"],
            "long_list": [
                "this",
                "is",
                "a",
                "rather",
                "long\nlist",
                "and should be broken up because of its width",
            ],
            "non_ascii": "汉语",
            1: 2,
        }
        json.dump(data, sys.stdout, cls=CompactJSONEncoder, ensure_ascii=False)
        exit()

    json.dump(json.load(sys.stdin), sys.stdout, cls=CompactJSONEncoder)
