from __future__ import annotations

from abc import abstractmethod

from .base import ResultStore


class CloudResultStore(ResultStore):
    @abstractmethod
    def save(self, result):
        ...
        