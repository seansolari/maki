
from abc import ABC, abstractmethod
from typing import Any, Type

from maki.benchmark.data.models import Dataset, DatasetReference


class DatasetMaterialiser(ABC):
    registry = {}
    
    def __init_subclass__(cls, **kwargs) -> None:
        super().__init_subclass__(**kwargs)

        dtype = getattr(cls, "input_type", None)
        
        if not dtype:
            raise TypeError(f"DatasetMaterialiser {cls} does not declare required `input_type` property.")
        elif not isinstance(dtype, Type):
            raise TypeError(f"DatasetMaterialiser {cls} `input_type` does not point to a type.")
        
        DatasetMaterialiser.registry[dtype] = cls
    
    @classmethod
    def from_type(cls, dtype: Type["DatasetReference"]) -> Type["DatasetMaterialiser"]:
        try:
            return DatasetMaterialiser.registry[dtype]
        except KeyError:
            raise TypeError(f"Unrecognosed dataset type: {dtype.name}")
    
    @classmethod
    def from_dataset(cls, dataset: "DatasetReference") -> Type["DatasetMaterialiser"]:
        return DatasetMaterialiser.from_type(type(dataset))
    
    @abstractmethod
    def materialise(self, token: Any) -> "Dataset":
        pass
