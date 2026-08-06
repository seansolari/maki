from __future__ import annotations

from pathlib import Path

from .cache import DatasetCache
from .generators import SyntheticGeneratorRegistry
from .materializer import (
    DatasetMaterializerRegistry,
    DerivedDatasetMaterializer,
    LocalDatasetMaterializer,
    RealDatasetMaterializer,
    SyntheticDatasetMaterializer,
)
from .base import DerivedDataset, LocalDataset, RealDataset, SyntheticDataset
from .registry import DatasetRegistry
from .transformations import TransformationRegistry


def create_default_dataset_infrastructure(
    cache_root: str | Path,
    dataset_registry: DatasetRegistry | None = None,
    generator_registry: SyntheticGeneratorRegistry | None = None,
    transformation_registry: TransformationRegistry | None = None,
) -> tuple[
    DatasetRegistry,
    DatasetMaterializerRegistry,
    DatasetCache,
    SyntheticGeneratorRegistry,
    TransformationRegistry,
]:
    dataset_registry = dataset_registry or DatasetRegistry()
    generator_registry = generator_registry or SyntheticGeneratorRegistry()
    transformation_registry = transformation_registry or TransformationRegistry()

    cache = DatasetCache(cache_root)
    materializers = DatasetMaterializerRegistry()

    materializers.register(
        LocalDataset,
        LocalDatasetMaterializer(cache),
    )

    materializers.register(
        RealDataset,
        RealDatasetMaterializer(cache),
    )

    materializers.register(
        SyntheticDataset,
        SyntheticDatasetMaterializer(
            cache=cache,
            generators=generator_registry,
        ),
    )

    materializers.register(
        DerivedDataset,
        DerivedDatasetMaterializer(
            cache=cache,
            dataset_registry=dataset_registry,
            materializer_registry=materializers,
            transformations=transformation_registry,
        ),
    )

    return (
        dataset_registry,
        materializers,
        cache,
        generator_registry,
        transformation_registry,
    )