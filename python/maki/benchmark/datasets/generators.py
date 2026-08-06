from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

from .errors import DatasetMaterializationError


class SyntheticDatasetGenerator(ABC):
    @property
    @abstractmethod
    def name(self) -> str:
        ...

    @abstractmethod
    def generate(
        self,
        output_path: Path,
        parameters: dict[str, Any],
        seed: int | None = None,
    ) -> Path:
        ...


class SyntheticGeneratorRegistry:
    def __init__(self) -> None:
        self._generators: dict[str, SyntheticDatasetGenerator] = {}

    def register(self, generator: SyntheticDatasetGenerator) -> None:
        if generator.name in self._generators:
            raise DatasetMaterializationError(
                f"Synthetic generator already registered: {generator.name}"
            )

        self._generators[generator.name] = generator

    def get(self, name: str) -> SyntheticDatasetGenerator:
        try:
            return self._generators[name]
        except KeyError as error:
            raise DatasetMaterializationError(
                f"Synthetic generator not registered: {name}"
            ) from error