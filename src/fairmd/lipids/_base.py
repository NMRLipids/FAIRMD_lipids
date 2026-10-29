"""Base abstract classes for the FairMD Lipid Databank."""

import copy
import sys
import typing
from abc import ABC, abstractmethod
from collections.abc import Iterable, MutableSet
from typing import Any, Generic, TypeVar

from tqdm import tqdm

from fairmd.lipids.bioschema_inheritance import resolve_block, template_stamp
from fairmd.lipids.molecules import Lipid, Molecule, NonLipid, lipids_set, solubles_set


def progress(
    iterable: Iterable | None = None,
    /,
    *,
    desc: str | None = None,
    total: int | None = None,
    disable: bool = False,
    **kwargs,
):
    """
    Wrap around tqdm with FAIRMD defaults.

    - Always writes to stdout
    - Progress is shown by default (including CI)
    - Can be explicitly disabled if needed
    """
    return tqdm(
        iterable or kwargs.pop("iter", None) or [],
        desc=desc,
        total=total,
        disable=disable,
        file=sys.stdout,
        **kwargs,
    )


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


_READ_ONLY_MSG = (
    "{kind} is a read-only view of a databank record's metadata. Assign through the record, "
    "e.g. system['bioschema_properties'] = block, or edit the copy that bioschema or "
    "to_readme() returns."
)


class FrozenDict(dict):
    """A ``dict`` that refuses changes; a read-only view handed out by :class:`MetadataProvider`.

    It is still a ``dict``: it compares equal to one and serialises to JSON as one.
    Copying it, shallow or deep, gives back a plain, mutable ``dict``.
    """

    def _read_only(self, *args, **kwargs):  # noqa: ANN002, ANN003, ANN202, ARG002
        raise TypeError(_READ_ONLY_MSG.format(kind=type(self).__name__))

    __setitem__ = __delitem__ = __ior__ = _read_only
    clear = pop = popitem = setdefault = update = _read_only

    def copy(self) -> dict:
        return dict(self)

    __copy__ = copy

    def __deepcopy__(self, memo: dict) -> dict:
        return {copy.deepcopy(k, memo): copy.deepcopy(v, memo) for k, v in self.items()}

    def __reduce__(self):  # noqa: ANN204
        return (dict, (dict(self),))


class FrozenList(list):
    """A ``list`` that refuses changes; see :class:`FrozenDict`."""

    def _read_only(self, *args, **kwargs):  # noqa: ANN002, ANN003, ANN202, ARG002
        raise TypeError(_READ_ONLY_MSG.format(kind=type(self).__name__))

    __setitem__ = __delitem__ = __iadd__ = __imul__ = _read_only
    append = extend = insert = pop = remove = clear = sort = reverse = _read_only

    def copy(self) -> list:
        return list(self)

    __copy__ = copy

    def __deepcopy__(self, memo: dict) -> list:
        return [copy.deepcopy(v, memo) for v in self]

    def __reduce__(self):  # noqa: ANN204
        return (list, (list(self),))


def freeze(value):  # noqa: ANN001, ANN201
    """A read-only copy of ``value``: dicts and lists, at every depth, become frozen."""
    if isinstance(value, dict):
        return FrozenDict({key: freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return FrozenList(freeze(item) for item in value)
    return value


class MetadataProvider(ABC):
    """Abstract base class for a databank record that carries metadata.

    A record's README stores a *compact* ``bioschema_properties`` block: what the
    ``bioschema_properties.yaml`` template of the record's kind already states is
    left out. Every public accessor returns the *resolved* block instead -- the
    kind's template with the record's own block applied on top (two layers; one
    template per kind, in :meth:`bioschema_template_dir`; merge rules in
    :mod:`fairmd.lipids.bioschema_inheritance`) -- so these are always equal::

        record.metadata["bioschema_properties"] == record.bioschema == record["bioschema_properties"]

    The views are read-only, so a change can never be made to a copy and lost:

    - ``record["bioschema_properties"]`` and the block inside :attr:`metadata`
      are frozen at every depth (:class:`FrozenDict`, :class:`FrozenList`).
      Change the block by assigning a new one to ``record["bioschema_properties"]``.
    - :attr:`metadata` itself is frozen at its top level. Its other values are
      the record's own objects, so changing one in place changes the record,
      exactly as ``record[key]`` does.
    - :attr:`bioschema` is a new, mutable ``dict``: a value to use, not a view.

    Only code inside the class sees the stored form, through
    :attr:`_stored_metadata`. The one sanctioned way out is :meth:`to_readme`,
    which is what any code writing a ``README.yaml`` back must use; a frozen view
    cannot be written by a YAML dumper, which keeps the template from being baked
    into a file.
    """

    BIOSCHEMA_KEY = "bioschema_properties"

    @property
    @abstractmethod
    def _stored_metadata(self) -> dict:
        """The README as loaded, compact block included. For use inside the class only."""

    @classmethod
    @abstractmethod
    def bioschema_template_dir(cls) -> str:
        """The folder holding the ``bioschema_properties.yaml`` every record of this kind inherits."""

    def _resolved_block(self) -> FrozenDict | None:
        """The resolved block, frozen; None when the README has none.

        Cached per record and recomputed only when the stored block is replaced
        or the template file changes, so reading the block -- or any key of
        :attr:`metadata` -- does not re-read and re-copy the template each time.
        """
        stored = self._stored_metadata.get(self.BIOSCHEMA_KEY)
        if stored is None:
            return None
        template_dir = self.bioschema_template_dir()
        stamp = (str(template_dir), template_stamp(template_dir))
        cached = self.__dict__.get("_bioschema_cache")
        if cached is not None and cached[0] is stored and cached[1] == stamp:
            return cached[2]
        resolved = freeze(resolve_block(stored, template_dir))
        self._bioschema_cache = (stored, stamp, resolved)
        return resolved

    def _metadata_item(self, key: str):  # noqa: ANN202
        """One metadata value as the public sees it: the block resolved, anything else as stored."""
        if key == self.BIOSCHEMA_KEY and self._stored_metadata[key] is not None:
            return self._resolved_block()
        return self._stored_metadata[key]

    @property
    def metadata(self) -> FrozenDict:
        """The record's metadata, with its Bioschemas block resolved; read-only at its top level.

        ``bioschema_properties`` is present only when the README has a block.
        """
        return FrozenDict({key: self._metadata_item(key) for key in self._stored_metadata})

    @property
    def bioschema(self) -> dict | None:
        """The record's complete Bioschemas ``Dataset`` block.

        :return: the kind's template updated by the record's own block, as a new,
                 mutable dict; None when the record carries no block yet.
        """
        resolved = self._resolved_block()
        return None if resolved is None else copy.deepcopy(resolved)

    def to_readme(self) -> dict:
        """The dict to write to ``README.yaml``: the stored form, with the compact block.

        Nothing inherited is in it, so writing it back never bakes the kind's
        template into the file. A deep copy: changing it does not change the record.
        """
        return copy.deepcopy(dict(self._stored_metadata))


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
