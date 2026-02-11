"""Base abstract classes for the FairMD Lipid Databank."""

import typing
from abc import ABC, abstractmethod
from collections.abc import MutableSet
from typing import Any, Generic, TypeVar

from fairmd.lipids.molecules import Lipid, Molecule, NonLipid, lipids_set, solubles_set


class SampleComposition(ABC):
    """Abstract base class representing a sample composition in the databank."""

    _content: dict[str, Molecule]

    def __init__(self) -> None:
        """Initialize the SampleComposition object."""
        self._content = {}

    @abstractmethod
    def _initialize_content(self) -> None:
        """Initialize the content of the sample composition."""

    @property
    def content(self) -> dict[str, Molecule]:
        """Returns dictionary of molecule objects."""
        return self._content

    @property
    def lipids(self) -> dict[str, Lipid]:
        """Returns dictionary of lipid molecule objects."""
        return {k: v for k, v in self.content.items() if k in lipids_set}

    @property
    def solubles(self) -> dict[str, NonLipid]:
        """Returns dictionary of non-lipid molecule objects."""
        return {k: v for k, v in self.content.items() if k in solubles_set and k != "SOL"}

    @abstractmethod
    def membrane_composition(self, basis: typing.Literal["molar", "mass"] = "molar") -> dict[str, float]:
        """Return the composition of the membrane in system.

        :param which: Type of composition to return. Options are:
                      - "molar": compute molar fraction
                      - "mass": compute mass fraction
        :return: dictionary (universal molecule name -> value)
        """

    @abstractmethod
    def get_hydration(self, basis: typing.Literal["number", "mass"] = "number") -> float:
        """Get system hydration."""

    @abstractmethod
    def solution_composition(self, basis: typing.Literal["molar", "mass"] = "molar") -> dict[str, float]:
        """Return the composition of the solution in system."""


T = TypeVar("T")


class CollectionSingleton(MutableSet[T], ABC, Generic[T]):
    """A generic, mutable set collection for databank items."""

    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        """Initialize the Empty Collection."""
        self._items: list[T] = list()
        self._ids: set = set()

    @classmethod
    def clear_instance(cls):
        """Clear the singleton instance. For testing purposes only."""
        cls._instance = None

    @abstractmethod
    def _test_item_type(self, item: Any) -> bool:
        """Test if an item is of the proper type for the collection."""

    @abstractmethod
    def _get_item_id(self, item: T) -> str:
        """Get the unique identifier of an item."""

    def __contains__(self, item: Any) -> bool:
        """Check if an item is in the set by instance or by ID."""
        if isinstance(item, str):
            item = item.upper()
        return (self._test_item_type(item) and item in self._items) or (item in self._ids)

    def __iter__(self):
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, index: int) -> T:
        return self._items[index]

    def add(self, item: T) -> None:
        """Add an item to the set."""
        if self._test_item_type(item):
            self._items.append(item)
            id = self._get_item_id(item)
            if isinstance(id, str):
                id = id.upper()
            if id in self._ids:
                msg = f"Item with ID '{id}' already exists in {type(self).__name__}."
                raise KeyError(msg)
            self._ids.add(id)
        else:
            msg = f"Only proper instances can be added to {type(self).__name__}."
            raise TypeError(msg)

    def discard(self, id: str | int) -> None:
        """Remove an item from the set without raising an error if it does not exist."""
        raise NotImplementedError("This method should be implemented for non-set.")
        # if isinstance(id, str):
        #     id = id.upper()
        # item_to_remove = self.get(id)
        # if item_to_remove:
        #     # self._items.discard(item_to_remove) MUST BE IMPLEMENTED FOR LIST
        #     self._ids.discard(id)

    def get(self, key: str | int, default: Any = None) -> T | None:
        """Get an item by its ID (case-insensitive)."""
        if isinstance(key, str):
            key = key.upper()
        if key in self._ids:
            for item in self._items:
                comparison_id = self._get_item_id(item)
                if isinstance(comparison_id, str):
                    comparison_id = comparison_id.upper()
                if comparison_id == key:
                    return item
        return default

    def __repr__(self) -> str:
        return f"{type(self).__name__}({sorted(list(self._ids))})"

    @property
    def ids(self) -> set:
        """The set of unique identifiers for all items in the collection."""
        return self._ids

    @staticmethod
    def load_from_data() -> "CollectionSingleton":
        """Load collection data from the designated directory."""
        msg = "This method should be implemented in subclasses."
        raise NotImplementedError(msg)

    # TODO: schedule for removing
    def loc(self, key: str | int) -> T:
        """Locate an item by its unique identifier."""
        return self.get(key)
