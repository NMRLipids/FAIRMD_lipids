"""Helpers for combining metadata supplied by databank objects."""

from abc import ABC, abstractmethod
from collections.abc import Mapping
from copy import deepcopy


class MetadataProvider(ABC):
    """Base class for objects exposing inherited and local metadata.

    Implementations only provide the two sources; the public mapping is always
    produced in one place and is therefore consistent for all providers.
    """

    def __init__(self) -> None:
        self._metadata: dict | None = None

    @property
    @abstractmethod
    def _parent_metadata(self) -> Mapping:
        """Metadata inherited from the provider's parent/context."""

    @property
    @abstractmethod
    def _current_metadata(self) -> Mapping:
        """Metadata belonging to this provider."""

    @property
    def metadata(self) -> Mapping:
        """Public metadata exposed by the provider."""
        if self._metadata is None:
            self._metadata = self._merge_mappings(self._parent_metadata, self._current_metadata)
        return self._metadata

    @staticmethod
    def _merge_mappings(parent: Mapping, child: Mapping) -> dict:
        merged = deepcopy(dict(parent))
        for key, value in child.items():
            if isinstance(value, Mapping) and isinstance(merged.get(key), Mapping):
                merged[key] = MetadataProvider._merge_mappings(merged[key], value)
            else:
                merged[key] = deepcopy(value)
        return merged
