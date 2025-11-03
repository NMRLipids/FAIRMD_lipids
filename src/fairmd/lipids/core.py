"""
Core databank class and system initialization function.

Imported by :mod:`fairmd.lipids.databankLibrary` by default.
Can be imported without additional libraries to scan Databank system file tree!
"""

import os
import sys
from collections.abc import Iterable, Iterator, Mapping, MutableMapping, Sequence

import yaml

from fairmd.lipids import FMDL_SIMU_PATH
from fairmd.lipids.settings.molecules import Lipid, Molecule, NonLipid, lipids_set, molecules_set


class System(MutableMapping):
    """
    Main Databank single object.

    It is an extension of a dictionary with additional functionality.
    """

    def __init__(self, data: dict | Mapping) -> None:
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
        elif isinstance(data, Mapping):
            self._store.update(dict(data))
        else:
            expect_type_msg = "Expected dict or Mapping"
            raise TypeError(expect_type_msg)

        self._content = {}
        for k, v in self["COMPOSITION"].items():
            mol = None
            if k in lipids_set:
                mol = Lipid(k)
            elif k in molecules_set:
                mol = NonLipid(k)
            else:
                mol_not_found_msg = f"Molecule {k} is not in the set of lipids or molecules."
                raise ValueError(mol_not_found_msg)
            mol.register_mapping(v["MAPPING"])
            self._content[k] = mol

    def __getitem__(self, key: str):  # noqa: ANN204
        return self._store[key]

    def __setitem__(self, key: str, value) -> None:
        self._store[key] = value

    def __delitem__(self, key: str) -> None:
        del self._store[key]

    def __iter__(self) -> Iterator:
        return iter(self._store)

    def __len__(self) -> int:
        return len(self._store)

    @property
    def readme(self) -> dict:
        """Get the README dictionary of the system in true dict format.

        :return: dict-type README (dict)
        """
        return self._store

    @property
    def content(self) -> dict[str, Molecule]:
        """Returns dictionary of molecule objects."""
        return self._content

    def __repr__(self) -> str:
        return f"System({self._store['ID']}): {self._store['path']}"


class SystemsCollection(Sequence[System]):
    """Immutable collection of system dicts. Can be accessed by ID using loc()."""

    def __init__(self, iterable: Iterable[System] = []) -> None:
        self._data = iterable
        self.__get_index_byid()

    def __get_index_byid(self) -> None:
        self._idx: dict = {}
        for i in range(len(self)):
            if "ID" in self[i]:
                self._idx[self[i]["ID"]] = i

    def __getitem__(self, i: int) -> System:
        return self._data[i]

    def __len__(self) -> int:
        return len(self._data)

    def loc(self, sid: int) -> System:
        """Locate system by its ID.

        :param sid: System ID
        :return: System object with ID `sid`
        """
        return self._data[self._idx[sid]]


class Databank:
    """
    Representation of all simulation in the NMR lipids databank.

    `path` should be the local location of `{FMDL_DATA_PATH}/Simulations/` in
    the FAIRMD Lipids folder. Example usage to loop over systems:

    .. code-block:: python

        path = 'BilayerData/Simulations/'
        db_data = databank(path)
        systems = db_data.get_systems()

        for system in systems:
            print(system)
    """

    def __init__(self) -> None:
        self.path = FMDL_SIMU_PATH
        __systems = self.__load_systems__()
        self._systems: SystemsCollection = SystemsCollection(__systems)
        print("Databank initialized from the folder:", os.path.realpath(self.path))

    def __load_systems__(self) -> list[System]:
        systems: list[System] = []
        rpath = os.path.realpath(self.path)
        for subdir, _dirs, files in os.walk(rpath):
            for filename in files:
                filepath = os.path.join(subdir, filename)
                if filename == "README.yaml":
                    ydict = {}
                    try:
                        with open(filepath) as yaml_file:
                            ydict.update(yaml.load(yaml_file, Loader=yaml.FullLoader))
                        content = System(ydict)
                    except (FileNotFoundError, PermissionError) as e:
                        sys.stderr.write(f"""
!!README LOAD ERROR!!
Problems while loading on of the files required for the system: {e}
System path: {subdir}
System: {ydict!s}\n""")
                    except Exception as e:
                        sys.stderr.write(f"""
!!README LOAD ERROR!!
Unexpected error: {e}
System: {ydict!s}\n""")
                    else:
                        relpath = os.path.relpath(filepath, rpath)
                        content["path"] = relpath[:-11]
                        systems.append(content)
        return systems

    def get_systems(self) -> SystemsCollection:
        """List all systems in the FAIRMD Lipids."""
        return self._systems


def initialize_databank() -> SystemsCollection:
    """
    Intialize the FAIRMD Lipids.

    :return: list of dictionaries that contain the content of README.yaml files for
             each system.
    """
    db_data = Databank()
    return db_data.get_systems()


# TODO: is not used at all in the project!!
def print_README(system: str | Mapping) -> None:  # noqa: N802
    """
    Print the content of ``system`` dictionary in human readable format.

    :param system: FAIRMD Lipids dictionary defining a simulation or "example".
    """
    if system == "example":
        current_folder = os.path.dirname(os.path.realpath(__file__))
        readme_path = os.path.join(current_folder, "settings", "READMEexplanations.yaml")
        with open(readme_path) as file:
            readme_file = yaml.safe_load(file)
    else:
        readme_file = system

    for key in readme_file:
        print("\033[1m" + key + ":" + "\033[0m")
        print(" ", readme_file[key])
