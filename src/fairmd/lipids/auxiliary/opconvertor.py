"""
Auxiliary functions for converting OP data.

Helps to build a nicely formatted OP dictionary from raw OP data.
Fragmentation is handled.
"""

import re

from natsort import natsorted

from fairmd.lipids.molecules import Lipid


class NamingRegistry:
    """Registry for naming conventions."""

    _registry: list = []

    @classmethod
    def _register(cls, name: str, func) -> None:
        """Register a naming convention function.

        :param name: Name of the fragment.
        :param func: Function implementing the convention.
        """
        cls._registry += [(name, func)]

    @classmethod
    def apply(cls, opdic: dict) -> None:
        """Make whatever required to apply naming conventions to fragmentized dictionary.

        :param opdic: Fragmented dictionary.
        """
        if not cls._registry:
            cls._initialize()
        cls._apply_naming(opdic)
        cls._apply_sorting(opdic)

    @classmethod
    def _apply_sorting(cls, opdic: dict) -> None:
        """Sort every fragment list by C atom number."""
        for frag_name in opdic:
            opdic[frag_name] = natsorted(opdic[frag_name], key=lambda x: x["C"] + "__" + x["H"])

    @classmethod
    def _apply_naming(cls, opdic: dict) -> None:
        """Apply naming conventions to the fragmented dictionary."""
        for frag_func in cls._registry:
            frag_name, func = frag_func
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
        # usual FA (glycerolipids)
        def _snX_c_renamer(row: dict) -> dict:  # noqa: N802
            match = re.match(r"M_G[12]C([0-9]{1,2})_M", row["C"])
            if match and len(match.groups()) == 1:
                idx = int(match[1])
                row["C"] = str(idx - 1)
            return row

        cls._register("sn-1", _snX_c_renamer)
        cls._register("sn-2", _snX_c_renamer)

        # usual FA (cardiolipins)
        def _snXX_c_renamer(row: dict) -> dict:  # noqa: N802
            match = re.match(r"M_G[12][12]C([0-9]{1,2})_M", row["C"])
            if match and len(match.groups()) == 1:
                idx = int(match[1])
                row["C"] = str(idx - 1)
            return row

        cls._register("sn-1 1", _snXX_c_renamer)
        cls._register("sn-1 2", _snXX_c_renamer)
        cls._register("sn-2 1", _snXX_c_renamer)
        cls._register("sn-2 2", _snXX_c_renamer)

        # FA (sphingolipids)
        def _fa_c_renamer(row: dict) -> dict:
            match = re.match(r"M_N1C([0-9]{1,2})_M", row["C"])
            if match and len(match.groups()) == 1:
                idx = int(match[1])
                row["C"] = str(idx)
            return row

        cls._register("fa", _fa_c_renamer)

        # glycerol backbone (glycerolipids)
        def _gbb_c_renamer(row: dict) -> dict:
            match = re.match(r"M_G([1-3])_M", row["C"])
            if match and len(match.groups()) == 1:
                idx = int(match[1])
                row["C"] = f"g{idx}"
            return row

        cls._register("glycerol backbone", _gbb_c_renamer)

        # glycerol bb (cardiolipins)
        def _gbbX_c_renamer(row: dict) -> dict:
            match = re.match(r"M_G[1-2]([1-3])_M", row["C"])
            if match and len(match.groups()) == 1:
                idx = int(match[1])
                row["C"] = f"g{idx}"
            return row

        cls._register("glycerol backbone 1", _gbbX_c_renamer)
        cls._register("glycerol backbone 2", _gbbX_c_renamer)

        # usual headgroup (phospholipids)
        def _headgroup_c_renamer(row: dict) -> dict:
            if row["C"] == "M_G3C4_M":
                row["C"] = "α"
            elif row["C"] == "M_G3C5_M":
                row["C"] = "β"
            elif re.match(r"M_G3N6C[1-3]_M", row["C"]):
                row["C"] = "γ"
            return row

        cls._register("headgroup", _headgroup_c_renamer)

        # intermediate glycerol (CL)
        def _midglyc_c_renamer(row: dict) -> dict:
            match = re.match(r"M_G0([1-3])_M", row["C"])
            if match and len(match.groups()) == 1:
                idx = int(match[1])
                row["C"] = f"g{idx}'"
            return row

        cls._register("headgroup", _midglyc_c_renamer)

        # universal H renamer
        def _h_renamer(row: dict) -> dict:
            match = re.match(r"M_.+H([1-4])_M", row["H"])
            if match and len(match.groups()) == 1:
                idx = int(match[1])
                row["H"] = str(idx)
            return row

        cls._register("_all_", _h_renamer)

        # universal plain C renamer
        def _plain_c_renamer(row: dict) -> dict:
            match = re.match(r"M_C([0-9]+)_M", row["C"])
            if match and len(match.groups()) == 1:
                idx = int(match[1])
                row["C"] = str(idx)
            return row

        cls._register("_all_", _plain_c_renamer)


def build_nice_OPdict(src: dict, lipid: Lipid) -> dict:
    """Build nicely formatted OP dictionary from raw OP data.

    Handles fragmentation of lipids.

    :param src: raw OP data (dict)
    :param lipid: Lipid object
    :return: nicely formatted OP dictionary
    """

    # Helper function to convert NaN to None for better
    # JSON compatibility in output
    def _rnan(x: float) -> float | None:
        return None if x != x else x

    def _fragmentize(src: dict, mdict: dict) -> dict:
        r = {}
        for apair, opvals in src.items():
            atom_c, atom_h = apair.split(" ")
            if atom_c not in mdict:
                msg = f"Atom {atom_c} not found in mapping dictionary."
                raise ValueError(msg)
            frag_c = mdict[atom_c].get("FRAGMENT", "total")
            if frag_c not in r:
                r[frag_c] = []
            r[frag_c].append(
                {
                    "C": atom_c,
                    "H": atom_h,
                    "OP": opvals[0],
                    "STD": _rnan(opvals[1]) if len(opvals) > 1 else None,
                },
            )
            r[frag_c].sort(key=lambda x: x["C"])
        return r

    nice_OPdict: dict = _fragmentize(src, lipid.mapping_dict)
    # rename C and H atoms for registry fragments
    NamingRegistry.apply(nice_OPdict)
    return nice_OPdict
