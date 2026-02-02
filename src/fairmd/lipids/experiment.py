"""Module for handling experimental data entries in the databank.

This module defines classes to represent different types of experimental datasets,
they inherit from a common abstract base class :class:`Experiment` which implements
:class:`SampleComposition <fairmd.lipids._base.SampleComposition>` for retrieving sample
composition information.

Two concrete experiment types are implemented: :class:`OPExperiment` for ssNMR order parameter data
and :class:`FFExperiment` for form factor data.
"""

import json
import os
import typing
from abc import abstractmethod
from typing import Any, Literal

import numpy as np
import yaml

from fairmd.lipids import FMDL_EXP_PATH
from fairmd.lipids._base import CollectionSingleton, SampleComposition
from fairmd.lipids.molecules import Lipid, NonLipid


class ExperimentError(Exception):
    """Experiment-related exception"""


class Experiment(SampleComposition):
    """Abstract base class representing an experimental dataset in the databank."""

    _exp_id: str
    _metadata: dict | None = None
    _data: dict | None = None

    def __init__(self, exp_id: str, path: str):
        """
        Initialize the Experiment object.

        :param exp_id: The unique identifier for the experiment.
        :param path: The absolute path to the experiment's directory.
        """
        self._exp_id = exp_id
        self._path = path
        self._populate_meta_data()
        self._initialize_content()

    def _get_path(self) -> str:
        """Return the absolute path to the experiment's directory."""
        return self._path

    def _populate_meta_data(self) -> None:
        """Populate metadata from the README.yaml file."""
        self._metadata = {}
        meta_path = os.path.join(self._get_path(), "README.yaml")
        if os.path.isfile(meta_path):
            with open(meta_path) as yaml_file:
                self._metadata = yaml.load(yaml_file, Loader=yaml.FullLoader)
            # remove None entries
            _keys_to_check = list(self._metadata.keys())
            for k in _keys_to_check:
                if self._metadata[k] is None:
                    del self._metadata[k]
        else:
            msg = f"Metadata file (README.yaml) not found for experiment '{self._exp_id}'."
            raise FileNotFoundError(msg)

    @property
    def metadata(self) -> dict:
        """Access the experiment's metadata."""
        if self._metadata is None:
            self._populate_meta_data()
        return self._metadata

    def __getitem__(self, key: str):
        return self.metadata[key]

    def __repr__(self):
        return f"{type(self).__name__}(id='{self.exp_id}')"

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.exp_id == other.exp_id

    def __hash__(self):
        return hash(self.exp_id)

    @property
    def readme(self) -> dict:
        """Provides access to the experiment's metadata (for backward compatibility)."""
        return self.metadata

    @property
    def exp_id(self) -> str:
        """The unique identifier of the experiment."""
        return self._exp_id

    @property
    def path(self) -> str:
        """The absolute path to the experiment's directory."""
        return self._path

    @property
    @abstractmethod
    def data(self) -> dict:
        """Provide access to the experiment's data."""

    @abstractmethod
    def verify_data(self) -> None:
        """Verify the integrity and consistency of the experiment's data."""

    @property
    @abstractmethod
    def exptype(self) -> str:
        """The type of the experiment."""

    @classmethod
    def target_folder(cls) -> str:
        """Get target folder name for the experiment type."""
        msg = "This method should be implemented in subclasses."
        raise NotImplementedError(msg)

    # Implementation of SampleComposition interface

    def _initialize_content(self) -> None:
        null_concentration_threshold = 1e-7  # 0.1 uM
        self._content = {}
        for k in self.metadata["MOLAR_FRACTIONS"]:
            lip = Lipid(k)
            lip.register_mapping()
            self._content[k] = lip
        for k, v in self.metadata.get("ION_CONCENTRATIONS", {}).items():
            if np.abs(float(v)) < null_concentration_threshold:
                continue
            self._content[k] = NonLipid(k)
        for k in self.metadata.get("COUNTER_IONS", {}):
            if k not in self._content:
                self._content[k] = NonLipid(k)

    def get_hydration(self, basis: typing.Literal["number", "mass"] = "number") -> float:
        if basis == "mass":
            return float(self.metadata["TOTAL_HYDRATION"])
        if basis == "number":
            tlc = self.metadata["TOTAL_LIPID_CONCENTRATION"]
            if tlc == "full hydration":
                return 70
            return 55.5 / float(tlc)  # water per lipid from outdated field
        msg = "Basis must be 'molar' or 'mass'"
        raise ValueError(msg)

    def membrane_composition(self, basis: typing.Literal["molar", "mass"] = "molar") -> dict[str, float]:
        if basis == "molar":
            return self.metadata["MOLAR_FRACTIONS"]
        if basis == "mass":
            comp: dict[str, float] = {}
            total_mass = 0.0
            for k, v in self.metadata["MOLAR_FRACTIONS"].items():
                mol: Lipid = self._content[k]
                mw = float(mol.metadata["bioschema_properties"]["molecularWeight"])
                comp[k] = v * mw
                total_mass += comp[k]
            for k in comp:
                comp[k] /= total_mass
            return comp
        msg = "Basis must be 'molar' or 'mass'"
        raise ValueError(msg)

    def solution_composition(self, basis="molar"):
        if basis == "molar":
            concs = {}
            _tlc = 55.5 / self.get_hydration("number")
            for k, v in self.metadata.get("ION_CONCENTRATIONS", {}).items():
                if k not in self.solubles:
                    continue
                concs[k] = float(v)
            for k, v in self.metadata.get("COUNTER_IONS", {}).items():
                if k not in concs:
                    concs[k] = 0.0
                concs[k] += _tlc * self.membrane_composition("molar")[v]
            return concs
        if basis == "mass":
            msg = "Mass basis for solution composition is not implemented yet"
            raise NotImplementedError(msg)
        msg = "Basis must be 'molar' or 'mass' (not implemented yet)"
        raise ValueError(msg)


class OPExperiment(Experiment):
    """Represents an ssNMR order parameter experiment.

    :attr:`data` provides access to a dictionary where keys are lipid names and values are dictionaries
    mapping atom unique names to their order parameters.

    :meth:`verify_data` checks that all atom unique names in the data correspond to known atoms in the lipid mapping.
    """

    @property
    def data(self) -> dict:
        if self._data is None:
            self._data = {}
            for fname in os.listdir(self._get_path()):
                if fname.endswith("_OrderParameters.json"):
                    molecule_name = fname.replace("_OrderParameters.json", "")
                    if molecule_name not in self.lipids:
                        msg = f"Data for non-existing molecule {molecule_name} in {self.exp_id}!"
                        raise ExperimentError(msg)
                    fpath = os.path.join(self._get_path(), fname)
                    with open(fpath) as json_file:
                        _tmpdic = json.load(json_file)
                        self._data[molecule_name] = {k: v[0] for k, v in _tmpdic.items()}
            if not self._data:
                msg = f"No order parameter data files found for experiment '{self.exp_id}'."
                raise ExperimentError(msg)
        return self._data

    def verify_data(self) -> None:
        for ln, opdic in self.data.items():
            lipid = self.lipids[ln]
            for uname_pair in opdic:
                if uname_pair.split(" ")[0] not in lipid.mapping_dict:
                    msg = (
                        "Order parameter data contains unknown atom "
                        f"'{uname_pair.split(' ')[0]}' for lipid '{ln}' in experiment '{self.exp_id}'."
                    )
                    raise ExperimentError(msg)

    @property
    def exptype(self) -> str:
        return "OrderParameters"

    @classmethod
    def target_folder(cls) -> str:
        return "OrderParameters"


class FFExperiment(Experiment):
    """Represents a form factor experiment."""

    @property
    def data(self) -> dict:
        if self._data is None:
            self._data = {}
            for fname in os.listdir(self._get_path()):
                if fname.endswith("_FormFactor.json"):
                    fpath = os.path.join(self._get_path(), fname)
                    with open(fpath) as json_file:
                        self._data = json.load(json_file)
                    break  # Assuming one form factor file per experiment
            if not self._data:
                msg = f"No form factor data file found for experiment '{self.exp_id}'."
                raise ExperimentError(msg)
        return self._data

    def verify_data(self) -> None:
        pass

    @property
    def exptype(self) -> str:
        return "FormFactors"

    @classmethod
    def target_folder(cls) -> str:
        return "FormFactors"


class ExperimentCollection(CollectionSingleton[Experiment]):
    """A collection of experiments."""

    def _test_item_type(self, item: Any) -> bool:
        return isinstance(item, Experiment)

    def _get_item_id(self, item: Experiment) -> str:
        return item.exp_id

    @staticmethod
    def load_from_data(exp_type: Literal["OPExperiment", "FFExperiment"] = "OPExperiment") -> "ExperimentCollection":
        """Load experiment data from the designated directory."""
        exp_types = {
            "OPExperiment": OPExperiment,
            "FFExperiment": FFExperiment,
        }
        if exp_type not in exp_types:
            msg = "..."
            raise ValueError(msg)
        collection = ExperimentCollection()
        for exp_cls in [exp_types[exp_type]]:
            path = os.path.join(FMDL_EXP_PATH, exp_cls.target_folder())
            if not os.path.isdir(path):
                continue
            for subdir, _, files in os.walk(path):
                if "README.yaml" in files:
                    # exp_id is the directory name relative to the exp_type directory
                    exp_id = os.path.relpath(subdir, path)
                    try:
                        exp = exp_cls(exp_id, subdir)
                        collection.add(exp)
                    except FileNotFoundError:
                        # Log this error?
                        pass
        return collection
