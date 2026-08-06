from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

from .errors import DatasetTransformationError


class DatasetTransformation(ABC):
    @property
    @abstractmethod
    def name(self) -> str:
        ...

    @abstractmethod
    def apply(
        self,
        input_path: Path,
        output_path: Path,
        parameters: dict[str, Any],
    ) -> Path:
        ...


class TransformationRegistry:
    def __init__(self) -> None:
        self._transformations: dict[str, DatasetTransformation] = {}

    def register(self, transformation: DatasetTransformation) -> None:
        if transformation.name in self._transformations:
            raise DatasetTransformationError(
                f"Transformation already registered: {transformation.name}"
            )

        self._transformations[transformation.name] = transformation

    def get(self, name: str) -> DatasetTransformation:
        try:
            return self._transformations[name]
        except KeyError as error:
            raise DatasetTransformationError(
                f"Transformation not registered: {name}"
            ) from error
