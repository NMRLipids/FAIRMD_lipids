"""
:module: settings/molecules.py

:description: Module file with definition of different global-level dictionaries.

There is a dictionary of lipids, ions, etc. If you add a lipid which is not yet
in the databank, you have to add it here!
"""

import os
import re
from abc import ABC, abstractmethod
from collections.abc import MutableSet
from typing import Any

import yaml

from fairmd.lipids import FMDL_MOL_PATH


class Molecule(ABC):
    """
    Abstract base class representing a molecule and its related operations.

    This class is designed to provide an interface for interacting with
    molecule-related files, which are stored in a molecule-related folder.
    It serves as a base for concrete implementations that need to define
    specific operations for handling molecule data.
    """

    _molname: str | None = None

    @abstractmethod
    def _get_path(self) -> str:
        """
        Return absolute path to molecule-related files.

        :return: str path
        """

    def register_mapping(self, fname: str) -> None:
        """
        Register mapping dictionary for the Molecule object

        :param fname: mapping filename (without path)
        :return:
        """
        self._mapping_fpath = os.path.join(self._get_path(), fname)
        assert os.path.isfile(self._mapping_fpath)

    @property
    def mapping_dict(self) -> dict:
        # preload on first call
        if self._mapping_dict is None:
            with open(self._mapping_fpath) as yaml_file:
                self._mapping_dict = yaml.load(yaml_file, Loader=yaml.FullLoader)
        return self._mapping_dict

    def __init__(self, name: str) -> None:
        """
        Create the molecule.

        :param name: The name of the molecule.
        :type name: str
        """
        self.__check_name(name)
        self._molname = name
        self._mapping_fpath = None
        self._mapping_dict = None

    @staticmethod
    def __check_name(name: str) -> None:
        """
        Check if the provided name contains only valid characters.

        This static method verifies that the input name string complies with
        the defined regular expression pattern, which restricts it to alphanumeric
        characters and underscores. If the name does not match the pattern, a
        ValueError is raised.

        :param name: A string representing the name to validate.

        :raises ValueError: If the name contains characters outside of the allowed set.
        """
        pat = r"[A-Za-z0-9_]+"
        if not re.match(pat, name):
            raise ValueError(f"Only {pat} symbols are allowed in Molecule name.")

    @abstractmethod
    def _populate_meta_data(self) -> None:
        """
        Populate metadata for the molecule from its YAML file.

        This method should be implemented in subclasses to load specific
        metadata related to the molecule.
        """

    @property
    @abstractmethod
    def metadata(self) -> dict:
        """Metadata for the molecule"""

    @property
    def name(self) -> str:
        """Molecule name"""
        return self._molname

    # comparison by name to behave in a set
    # It's case-insesitive as folder structure should work on mac/win

    def __eq__(self, other) -> bool:
        return isinstance(other, type(self)) and self.name.upper() == other.name.upper()

    def __hash__(self):
        return hash(self.name.upper())

    def __repr__(self):
        return f"{type(self).__name__}({self.name})"


class Lipid(Molecule):
    """
    Lipid class inherited from Molecule base. Contains all the molecules
    which belongs to the bilayer.
    """

    _lipids_dir: str = os.path.join(FMDL_MOL_PATH, "membrane")
    """Directory containing all lipid-subdirs"""

    def _get_path(self) -> str:
        return os.path.join(self._lipids_dir, self.name)

    def _populate_meta_data(self) -> None:
        """Populate metadata for the lipid from its YAML file."""
        self._metadata = {}
        meta_path = os.path.join(self._get_path(), "metadata.yaml")
        if os.path.isfile(meta_path):
            with open(meta_path) as yaml_file:
                self._metadata = yaml.load(yaml_file, Loader=yaml.FullLoader)
        else:
            raise FileNotFoundError(f"Metadata file not found for {self.name}.")

    @property
    def metadata(self) -> dict:
        """
        Return metadata for the lipid.

        :return: dict with metadata
        """
        if not hasattr(self, "_metadata"):
            self._populate_meta_data()
        return self._metadata

    def __init__(self, name: str) -> None:
        """
        Create the lipid.

        :param name: The name of the lipid.
        :type name: str
        """
        super().__init__(name)
        self._populate_meta_data()


class NonLipid(Molecule):
    """Class for non-bilayer molecules: solvent, ions, etc."""

    _nonlipids_dir: str = os.path.join(FMDL_MOL_PATH, "solution")
    """Directory containing all lipid-subdirs"""

    def _get_path(self) -> str:
        return os.path.join(self._nonlipids_dir, self.name)

    def _populate_meta_data(self) -> None:
        pass

    @property
    def metadata(self) -> dict:
        return {}


class MoleculeSet(MutableSet[Molecule], ABC):
    """Set repeating normal set functionality with some additional molecule-specific things."""

    def __init__(self, *args):
        """Initialize the LipidSet with optional initial elements."""
        self._items = set()
        self._names = set()
        for arg in args:
            self.add(arg)

    @abstractmethod
    def _test_item_type(self, item: Any) -> bool:
        """
        Test if item is of proper Molecule type

        :param item: any type variable
        :return: bool answer
        """

    @abstractmethod
    def _create_item(self, name: str) -> Molecule:
        """
        Construct the molecule of the proper type

        :param name: unique name of the molecule
        :return: constructed Molecule
        """

    def __contains__(self, item: Molecule):
        """Check if a lipid is in the set."""
        return (self._test_item_type(item) and item in self._items) or (
            isinstance(item, str) and item.upper() in self._names
        )

    def __iter__(self):
        return iter(self._items)

    def __len__(self):
        return len(self._items)

    def add(self, item: Molecule):
        """
        Add a lipid to the set.

        :param item: Can add either Molecule or str (then Molecule constructor
                     will be called)
        """
        if self._test_item_type(item):
            self._items.add(item)
            self._names.add(item.name.upper())
        elif isinstance(item, str):
            self._items.add(self._create_item(item))  # here we call Lipid constructor
            self._names.add(item.upper())
        else:
            raise TypeError(f"Only proper instances can be added to {type(self).__name__}.")

    def discard(self, item):
        """
        Remove a lipid from the set without raising an error if it does not exist.
        """
        if self._test_item_type(item):
            self._items.discard(item)
            self._names.discard(item.name.toupper())
        elif isinstance(item, str):
            if item.upper() not in self._names:
                return
            ifound = None
            for i in self._items:
                if i.name.upper() == item.upper():
                    ifound = i
                    break
            assert ifound is not None
            self._items.discard(ifound)
            self._names.discard(item.upper())

    def get(self, key: str, default=None) -> Molecule | None:
        """
        Get a molecule by its name.
        :param key: The name of the molecule to retrieve.
        :param default: The value to return if the molecule is not found.
        """
        if key.upper() in self._names:
            for item in self._items:
                if item.name.upper() == key.upper():
                    return item
        return default

    def __repr__(self):
        return f"{type(self).__name__}[{self._names}]"

    @property
    def names(self) -> set[str]:
        return self._names


class LipidSet(MoleculeSet):
    """MoleculeSet specialization for Lipid."""

    def _test_item_type(self, item: Any) -> bool:
        return isinstance(item, Lipid)

    def _create_item(self, name):
        return Lipid(name)

    @staticmethod
    def load_from_data():
        """
        Loads lipid data from the designated directory and returns a set of lipids.

        :rtype: LipidSet
        :return: An instance of loaded `LipidSet`.
        """
        lset = LipidSet()
        path = Lipid._lipids_dir
        sub_dirs = [name for name in os.listdir(path) if os.path.isdir(os.path.join(path, name))]
        for name in sub_dirs:
            lset.add(name)
        return lset


class NonLipidSet(MoleculeSet):
    """
    MoleculeSet specialization for NonLipid.
    """

    def _test_item_type(self, item: Any) -> bool:
        return isinstance(item, NonLipid)

    def _create_item(self, name):
        return NonLipid(name)

    @staticmethod
    def load_from_data():
        """
        Loads Nonlipid data from the designated directory and returns a set of lipids.

        :rtype: NonLipidSet
        :return: An instance of loaded `NonLipidSet`.
        """
        lset = NonLipidSet()
        path = NonLipid._nonlipids_dir
        sub_dirs = [name for name in os.listdir(path) if os.path.isdir(os.path.join(path, name))]
        for name in sub_dirs:
            lset.add(name)
        return lset


lipids_set: LipidSet = LipidSet.load_from_data()
""" MutableSet of possible lipids """

lipids_dict = lipids_set
""" @deprecated: Use lipids_set instead. """

molecules_set: NonLipidSet = NonLipidSet.load_from_data()
""" Dictionary of other than lipid molecules. """

molecules_dict = molecules_set
""" @deprecated: Use molecules_set instead"""

molecule_ff_set = {("FF" + x) for x in (lipids_set.names | molecules_set.names)}
"""
Dictionary containing possible force-field labels for molecules given by the contributor
 (used for README/info fields validation)
"""
