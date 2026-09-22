"""
Core databank class and system initialization function.

Imported by :mod:`fairmd.lipids.databankLibrary` by default.
Can be imported without additional libraries to scan Databank system file tree!
"""

import os
import sys
import typing
from collections.abc import MutableMapping
from typing import Any

import yaml

from fairmd.lipids import FMDL_SIMU_PATH
from fairmd.lipids._base import CollectionSingleton, SampleComposition
from fairmd.lipids.molecules import Lipid, NonLipid, lipids_set, solubles_set


class System(MutableMapping, SampleComposition):
    """
    Main Databank single object.

    It is an extension of a dictionary with additional functionality.
    """

    def __init__(self, data: dict | typing.Mapping) -> None:
        """
        Initialize the container for storing simulation record.

        :param data: README-type dictionary.
        :raises TypeError: If `data` is neither a `dict` nor another mapping type.
        :raises ValueError: If a molecule key in the "COMPOSITION" data does not
                            belong to the predefined set of lipids or molecules.
        """
        self._store: dict = {}
        if isinstance(data, dict):
            self._store.update(data)
        elif isinstance(data, typing.Mapping):
            self._store.update(dict(data))
        else:
            expect_type_msg = "Expected dict or Mapping"
            raise TypeError(expect_type_msg)

        self._initialize_content()

    def __getitem__(self, key: str):  # noqa: ANN204
        return self._store[key]

    def __setitem__(self, key: str, value) -> None:
        self._store[key] = value

    def __delitem__(self, key: str) -> None:
        del self._store[key]

    def __iter__(self) -> typing.Iterator:
        return iter(self._store)

    def __len__(self) -> int:
        return len(self._store)

    def __repr__(self) -> str:
        return f"System({self._store['ID']}): {self._store['path']}"

    @property
    def readme(self) -> dict:
        """Get the README dictionary of the system in true dict format.

        :return: dict-type README (dict)
        """
        return self._store

    @property
    def n_lipids(self) -> int:
        """Returns total number of lipid molecules in the system."""
        total = 0
        for k, v in self["COMPOSITION"].items():
            if k in lipids_set:
                total += sum(v["COUNT"])
        return total

    # Implementation of SampleComposition interface

    def _initialize_content(self) -> None:
        self._content = {}
        for k, v in self["COMPOSITION"].items():
            mol = None
            if k in lipids_set:
                mol = Lipid(k)
            elif k in solubles_set:
                mol = NonLipid(k)
            else:
                mol_not_found_msg = f"Molecule {k} is not in the set of lipids or molecules."
                raise ValueError(mol_not_found_msg)
            mol.register_mapping(v["MAPPING"])
            self._content[k] = mol

    def membrane_composition(self, basis: typing.Literal["molar", "mass"] = "molar") -> dict[str, float]:
        if basis not in ["molar", "mass"]:
            msg = "Basis must be 'molar' or 'mass'"
            raise ValueError(msg)
        comp: dict[str, float] = {}
        for k, v in self["COMPOSITION"].items():
            if k not in lipids_set:
                continue
            count = sum(v["COUNT"])
            comp[k] = count
        n_lipids = self.n_lipids
        if basis == "molar":
            for k in comp:
                comp[k] /= n_lipids
        else:  # (basis == "mass")
            total_mass = 0.0
            for k in comp:  # noqa: PLC0206 (modify dict while iterating)
                mol: Lipid = self._content[k]
                mw = mol.metadata["bioschema_properties"]["molecularWeight"]
                comp[k] *= mw
                total_mass += comp[k]
            for k in comp:
                comp[k] /= total_mass
        return comp

    def get_hydration(self, basis: typing.Literal["number", "mass"] = "number") -> float:
        if basis not in ["number", "mass"]:
            msg = "Basis must be 'molar' or 'mass'"
            raise ValueError(msg)
        if basis == "number":
            if "SOL" not in self["COMPOSITION"]:
                msg = "Cannot compute hydration for implicit water (system #{}).".format(self["ID"])
                raise ValueError(msg)
            hyval = self["COMPOSITION"]["SOL"]["COUNT"] / self.n_lipids
        else:  # basis == "mass"
            msg = "Mass hydration is not implemented yet."
            raise NotImplementedError(msg)
        return hyval

    def solution_composition(self, basis="molar"):
        if basis not in ["molar", "mass"]:
            msg = "Basis must be 'molar' or 'mass'"
            raise ValueError(msg)
        if not self.solubles:
            return {}  # pure water is allowed here (even for implicit water)
        if self["COMPOSITION"].get("SOL") is None:
            msg = "Cannot compute solution composition for implicit water (system #{}).".format(self["ID"])
            raise ValueError(msg)
        n_water = self["COMPOSITION"].get("SOL")
        comp: dict[str, float] = {}
        for k, v in self["COMPOSITION"].items():
            if k in lipids_set or k == "SOL":
                continue
            # convert to molar concentration
            # NOTE: we assume dilute solution. Sometimes not true!
            comp[k] = v["COUNT"] / n_water["COUNT"] * 55.5
        if basis == "molar":
            return comp
        # TODO: solubles doesn't have mass data yet
        msg = "Mass basis not implemented for solubles."
        raise NotImplementedError(msg)


class SystemsCollection(CollectionSingleton[System]):
    """Immutable collection of system dicts. Can be accessed by ID using loc()."""

    def _get_item_id(self, item: System) -> int:
        return item["ID"]

    def _test_item_type(self, item: Any) -> bool:
        return isinstance(item, System)

    @staticmethod
    def load_from_data() -> "SystemsCollection":
        """Load systems data from the designated directory."""
        print("Simulations are initialized from the folder:", os.path.realpath(FMDL_SIMU_PATH))
        systems = SystemsCollection()
        for subdir, _dirs, files in os.walk(FMDL_SIMU_PATH):
            for filename in files:
                filepath = os.path.join(subdir, filename)
                if filename == "README.yaml":
                    ydict = {}
                    try:
                        with open(filepath) as yaml_file:
                            ydict.update(yaml.load(yaml_file, Loader=yaml.FullLoader))
                    except (FileNotFoundError, PermissionError) as e:
                        sys.stderr.write(f"""
!!README LOAD ERROR!!
Problems while loading on of the files required for the system: {e}
System path: {subdir}
System: {ydict!s}\n""")
                    try:
                        content = System(ydict)
                    except Exception as e:
                        sys.stderr.write(f"""
!!README LOAD ERROR!!
Unexpected error: {e}
System: {ydict!s}\n""")
                    else:
                        content["path"] = os.path.relpath(subdir, FMDL_SIMU_PATH)
                        systems.add(content)
        return systems


def initialize_databank() -> SystemsCollection:
    """
    Returns Simulation collection (an alias).

    :return: list of dictionaries that contain the content of README.yaml files for
             each system.
    """
    return SystemsCollection.load_from_data()


# TODO: is not used at all in the project!!
def print_README(system: str | typing.Mapping) -> None:  # noqa: N802
    """
    Print the content of ``system`` dictionary in human readable format.

    :param system: FAIRMD Lipids dictionary defining a simulation or "example".
    """
    if system == "example":
        current_folder = os.path.dirname(os.path.realpath(__file__))
        readme_path = os.path.join(current_folder, "schema_validation", "schema", "READMEexplanations.yaml")
        with open(readme_path) as file:
            readme_file = yaml.safe_load(file)
    else:
        readme_file = system

    for key in readme_file:
        print("\033[1m" + key + ":" + "\033[0m")
        print(" ", readme_file[key])
