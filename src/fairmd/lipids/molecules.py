"""
:module: settings/molecules.py

:description: Module file with definition of different global-level dictionaries.

There is a dictionary of lipids, ions, etc. If you add a lipid which is not yet
in the databank, you have to add it here!
"""

import fnmatch
import os
import re
from abc import ABC, abstractmethod
from collections.abc import MutableSet
from typing import Any

import MDAnalysis as mda
import yaml

from fairmd.lipids import _HAS_RDKIT, FMDL_MOL_PATH


class MoleculeError(Exception):
    """
    Custom exception for molecule-related errors.

    @param message: Error message describing the issue.
    @param mol: Optional Molecule object related to the error.
    """

    def __init__(self, message: str, mol=None) -> None:
        msg = f"From molecule {mol}: {message}" if mol is not None else message
        self._mol = mol
        super().__init__(msg)


class MoleculeMappingError(MoleculeError):
    """Custom exception for molecule mapping errors."""

    def __init__(self, message: str, mol=None) -> None:
        if mol is None:
            msg = message
        elif mol.mapping_file is None:
            msg = f"From {mol}: {message}"
        else:
            msg = f"From {mol}[{mol.mapping_file}]: {message}" if mol is not None else message
        super().__init__(msg, mol=mol)


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
        self._disp_mapping = fname
        self._mapping_fpath = os.path.join(self._get_path(), fname)
        if not os.path.isfile(self._mapping_fpath):
            msg = f"Cannot find '{self._mapping_fpath}' mapping for molecule {self.name}"
            raise FileNotFoundError(msg)

    def check_mapping(self, u: mda.Universe, name: str) -> bool:
        """Check consistency of mapping file against Universe.

        :param u: MDAnalysis Universe
        :param name: string standing for residue name if it is not in the mapping file

        :raises MoleculeMappingError: if mapping is inconsistent
        """
        res_set = {v["RESIDUE"] for _, v in self.mapping_dict.items() if v.get("RESIDUE") is not None}
        res_set.add(name)
        res_atoms = u.select_atoms("resname " + " ".join(res_set))
        for a in res_atoms:
            _ = self.md2uan(a.name)
        for unm in self.mapping_dict:
            sel_str = self.uan2selection(unm, name)
            ats: mda.AtomGroup = u.select_atoms(sel_str)
            if ats.n_atoms == 0:
                msg = f"Atom {unm} was not found in the universe using selection '{sel_str}'."
                if name in sel_str:
                    msg += f" Could be mapping error or incorrect residue name '{name}'."
                raise MoleculeMappingError(msg, mol=self)

    @property
    def mapping_dict(self) -> dict:
        """Return mapping dictionary (load on first call)"""
        if self._mapping_fpath is None:
            msg = "Mapping file is not registered!"
            raise MoleculeError(msg, mol=self)
        if self._mapping_dict is None:
            try:
                with open(self._mapping_fpath) as yaml_file:
                    self._mapping_dict = yaml.safe_load(yaml_file)  # yaml.load(yaml_file, Loader=yaml.FullLoader)
            except OSError as e:
                msg = "Error opening mapping-file!"
                raise MoleculeError(msg, mol=self) from e
        return self._mapping_dict

    def md2uan(self, mdatomname: str, mdresname: str | None = None) -> str:
        """
        Convert MD atom name to the Universal Atom Name.

        :raises MoleculeMappingError: if the md atom name is not found in the mapping.
        :return: Universal Atom Name (str)
        """
        for universal_name, mrecord in self.mapping_dict.items():
            mapping_aname = mrecord["ATOMNAME"]
            # MDAnalysis uses fnmatch patterns for selection language
            # https://userguide.mdanalysis.org/stable/selections.html
            if "RESIDUE" in mrecord and mdresname is not None:
                mapping_rname = mrecord["RESIDUE"]
                if not fnmatch.fnmatch(mdresname, mapping_rname):
                    continue
            if fnmatch.fnmatch(mdatomname, mapping_aname):
                return universal_name
        emsg = f"Atom {mdatomname} is not found."
        raise MoleculeMappingError(emsg, mol=self)

    def uan2selection(self, uname: str, resname: str) -> str:
        """
        Convert the Universal Atom Name to MD atom name.

        :raises KeyError: if the universal name is not found in the mapping.
        :return: selection string for MDAnalysis
        """
        anm = self.mapping_dict[uname]["ATOMNAME"]
        selstr = f"name {anm}"
        if "RESIDUE" in self.mapping_dict[uname]:
            selstr += " and resname " + self.mapping_dict[uname]["RESIDUE"]
        else:
            selstr += " and resname " + resname
        return selstr

    def __init__(self, name: str) -> None:
        """
        Create the molecule.

        :param name: The name of the molecule.
        :type name: str
        """
        self.__check_name(name)
        self._molname = name
        self._disp_mapping = None
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

        :raises MoleculeError: If the name contains characters outside of the allowed set.
        """
        pat = r"[A-Za-z0-9_]+"
        if not re.match(pat, name):
            msg = f"Only {pat} symbols are allowed in Molecule name."
            raise MoleculeError(msg)

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

    @property
    def mapping_file(self) -> str:
        """Mapping file name"""
        return self._disp_mapping

    # comparison by name to behave in a set
    # It's case-insesitive as folder structure should work on mac/win

    def __eq__(self, other) -> bool:
        return isinstance(other, type(self)) and self.name.upper() == other.name.upper()

    def __hash__(self) -> int:
        return hash(self.name.upper())

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.name})"


class Lipid(Molecule):
    """
    Lipid class inherited from :class:`Molecule` base.

    Is used for all the molecules located in the bilayer.
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
            msg = f"Metadata file not found for {self.name}."
            raise FileNotFoundError(msg)

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

    def smi2uan(self, smile_id: int) -> str:
        """Convert SMILEIDX to universal atomname"""
        if smile_id < 0:
            msg = "ID<0: out of range"
            raise KeyError(msg)
        max_smid = 0
        for universal_name, mrecord in self.mapping_dict.items():
            smi = mrecord.get("SMILEIDX", -1)
            if smi == smile_id:
                return universal_name
            max_smid = max(max_smid, smi)
        if smile_id > max_smid:
            msg = f"ID>{max_smid}: out of range"
            raise KeyError(msg)

        emsg = f"Atom with SMILEIDX {smile_id} is not found."
        raise MoleculeMappingError(emsg, mol=self)

    @property
    def rdkit_object(self) -> "rdkit.Chem.Mol":  # noqa: F821
        """
        Return RDKit molecule object for the lipid.

        :return: RDKit molecule object
        """
        if not _HAS_RDKIT:
            msg = "RDKit is required for RDKit molecule object. Please install the 'rdkit' extra."
            raise ImportError(msg)

        from rdkit import Chem  # noqa: PLC0415

        if "smiles" not in self.metadata.get("bioschema_properties", {}):
            msg = (
                "SMILES is obligatory in `bioschema_properties->smiles` "
                "to request SMARTS and use cheminformatic integration."
            )
            raise MoleculeError(msg, mol=self)
        _smiles = self.metadata["bioschema_properties"]["smiles"]
        rdkit_mol = Chem.MolFromSmiles(_smiles)
        if rdkit_mol is None:
            msg = f"Invalid SMILES for molecule {self.name}: {_smiles}"
            raise MoleculeError(msg, mol=self)
        return rdkit_mol

    def atoms_by(self, query: str, id: int) -> list[str]:
        """
        Return list of universal atom names matching the SMARTS query|id.

        :param query: SMARTS pattern
        :type query: str
        :param id: Position inside SMARTS pattern
        :type id: int
        :return: list of universal atom names matching the query
        :rtype: list[str]
        """
        list_unames = []
        rdkit_mol = self.rdkit_object  # import check happens here

        from rdkit import Chem  # noqa: PLC0415

        patt = Chem.MolFromSmarts(query)
        if patt is None:
            msg = "SMART syntax error (see above)"
            raise ValueError(msg)
        if id < 0 or id > patt.GetNumHeavyAtoms() - 1:
            msg = f"SMARTS '{query}' has {patt.GetNumHeavyAtoms()} heavy atoms. {id} out of range."
            raise KeyError(msg)

        matches = rdkit_mol.GetSubstructMatches(patt)
        for match in matches:
            atom_idx = match[id]
            # USE FURHTER: rdk_atom = rdkit_mol.GetAtomWithIdx(atom_idx)
            list_unames.append(self.smi2uan(atom_idx))

        return list_unames


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

    def __contains__(self, item: Molecule) -> bool:
        """Check if a lipid is in the set."""
        return (self._test_item_type(item) and item in self._items) or (
            isinstance(item, str) and item.upper() in self._names
        )

    def __iter__(self):
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)

    def add(self, item) -> None:
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
            msg = f"Only proper instances can be added to {type(self).__name__}."
            raise TypeError(msg)

    def discard(self, item) -> None:
        """Remove a lipid from the set without raising an error if it does not exist."""
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
            if ifound is None:
                msg = f"Item '{item}' not found in the set."
                raise MoleculeError(msg)
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

    def __repr__(self) -> str:
        return f"{type(self).__name__}[{self._names}]"

    @property
    def names(self) -> set[str]:
        """Set of molecule names in the set."""
        return self._names


class LipidSet(MoleculeSet):
    """MoleculeSet specialization for Lipid."""

    def _test_item_type(self, item: Any) -> bool:
        return isinstance(item, Lipid)

    def _create_item(self, name: str) -> Molecule:
        return Lipid(name)

    @staticmethod
    def load_from_data() -> "LipidSet":
        """
        Load lipid data from the designated directory and returns a set of lipids.

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
    """MoleculeSet specialization for NonLipid."""

    def _test_item_type(self, item: Any) -> bool:
        return isinstance(item, NonLipid)

    def _create_item(self, name: str) -> Molecule:
        return NonLipid(name)

    @staticmethod
    def load_from_data() -> "NonLipidSet":
        """
        Load Nonlipid data from the designated directory and returns a set of lipids.

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
