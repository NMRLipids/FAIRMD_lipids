"""
Auxiliary functions for converting OP data.

Helps to build a nicely formatted OP dictionary from raw OP data.
Fragmentation is handled.
"""

import re

from fairmd.lipids.molecules import Lipid


class NamingRegistry:
    """Registry for naming conventions."""

    _registry: dict = {}

    @classmethod
    def _register(cls, name: str, func) -> None:
        """Register a naming convention function.

        :param name: Name of the fragment.
        :param func: Function implementing the convention.
        """
        cls._registry[name] = func

    @classmethod
    def apply(cls, opdic: dict):
        """Apply a naming convention functions.

        :param opdic: Fragmented dictionary.
        :return: Function implementing the convention.
        """
        if not cls._registry:
            cls._initialize()
        for frag_name, func in cls._registry.items():
            if frag_name in opdic:
                for i in range(len(opdic[frag_name])):
                    opdic[frag_name][i] = func(opdic[frag_name][i])
            elif frag_name == "_all_":
                for f in opdic:
                    for i in range(len(opdic[f])):
                        opdic[f][i] = func(opdic[f][i])

    # initialize the registry
    @classmethod
    def _initialize(cls):
        def _snX_c_renamer(row: dict) -> dict:
            match = re.match(r"M_G[12]C([0-9]{1,2})_M", row["C"])
            idx = int(match[1])
            row["C"] = str(idx - 1)
            return row

        cls._register("sn-1", _snX_c_renamer)

        def _gbb_c_renamer(row: dict) -> dict:
            match = re.match(r"M_G([1-3])_M", row["C"])
            idx = int(match[1])
            row["C"] = f"g{idx}"
            return row

        cls._register("glycerol backbone", _gbb_c_renamer)

        def _h_renamer(row: dict) -> dict:
            match = re.match(r"M_.+H([1-4])", row["H"])
            idx = int(match[1])
            row["H"] = str(idx)
            return row

        cls._register("_all_", _h_renamer)


def build_nice_OPdict(src: dict, lipid: Lipid) -> dict:
    """Build nicely formatted OP dictionary from raw OP data.

    Handles fragmentation of lipids.

    :param src: raw OP data (dict)
    :param lipid: Lipid object
    :return: nicely formatted OP dictionary
    """

    def _fragmentize(src, mdict):
        r = {}
        for apair, opvals in src.items():
            atom_c, atom_h = apair.split(" ")
            frag_c = mdict[atom_c].get("FRAGMENT", "total")
            if frag_c not in r:
                r[frag_c] = []
            r[frag_c].append(
                {"C": atom_c, "H": atom_h, "OP": opvals[0], "STD": opvals[1]},
            )
            r[frag_c].sort(key=lambda x: x["C"])
        return r

    nice_OPdict: dict = _fragmentize(src, lipid.mapping_dict)
    # rename C and H atoms for registry fragments
    NamingRegistry.apply(nice_OPdict)
    return nice_OPdict
