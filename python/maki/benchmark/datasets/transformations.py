from __future__ import annotations

from abc import ABC
from abc import abstractmethod
from typing import Any


class DatasetTransformation(ABC):

    @property
    @abstractmethod
    def name(self) -> str:
        ...

    @abstractmethod
    def apply(
        self,
        data: Any,
    ) -> Any:
        ...


class ReadSubsampleTransformation(DatasetTransformation):
    raise NotImplementedError()
