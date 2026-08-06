from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from urllib.request import urlretrieve

from .base import (
    BenchmarkDataset,
    DerivedDataset,
    LocalDataset,
    RealDataset,
    SyntheticDataset,
)
from .cache import DatasetCache
from .errors import (
    DatasetChecksumError,
    DatasetMaterializationError,
)
from .generators import SyntheticGeneratorRegistry
from .registry import DatasetRegistry
from .schema import MaterializedDataset
from .transformations import TransformationRegistry
from .utils import compute_checksum, copy_path


class DatasetMaterializer(ABC):
    @abstractmethod
    def materialize(
        self,
        dataset: BenchmarkDataset,
    ) -> MaterializedDataset:
        ...


class BaseDatasetMaterializer(DatasetMaterializer):
    def __init__(self, cache: DatasetCache) -> None:
        self.cache = cache

    def _materialized(
        self,
        dataset: BenchmarkDataset,
        path: Path,
    ) -> MaterializedDataset:
        return MaterializedDataset(
            dataset_id=dataset.dataset_id,
            version=dataset.version,
            path=str(path),
            metadata=dataset.metadata,
        )

    def _validate_checksum(
        self,
        dataset: BenchmarkDataset,
        path: Path,
    ) -> None:
        expected = dataset.metadata.checksum

        if not expected:
            return

        actual = compute_checksum(
            path,
            dataset.metadata.checksum_algorithm,
        )

        if actual != expected:
            raise DatasetChecksumError(
                f"Checksum mismatch for dataset {dataset.dataset_id}: "
                f"expected {expected}, got {actual}"
            )


class LocalDatasetMaterializer(BaseDatasetMaterializer):
    def materialize(
        self,
        dataset: BenchmarkDataset,
    ) -> MaterializedDataset:
        if not isinstance(dataset, LocalDataset):
            raise DatasetMaterializationError(
                "LocalDatasetMaterializer requires LocalDataset"
            )

        source = Path(dataset.path)

        if not source.exists():
            raise DatasetMaterializationError(
                f"Local dataset path does not exist: {source}"
            )

        self._validate_checksum(dataset, source)

        return self._materialized(dataset, source)


class RealDatasetMaterializer(BaseDatasetMaterializer):
    def materialize(
        self,
        dataset: BenchmarkDataset,
    ) -> MaterializedDataset:
        if not isinstance(dataset, RealDataset):
            raise DatasetMaterializationError(
                "RealDatasetMaterializer requires RealDataset"
            )

        cached_data_path = self.cache.data_path_for(dataset)

        if self.cache.exists(dataset):
            self._validate_checksum(dataset, cached_data_path)
            return self._materialized(dataset, cached_data_path)

        if not dataset.urls:
            raise DatasetMaterializationError(
                f"Real dataset has no URLs: {dataset.dataset_id}"
            )

        self.cache.prepare(dataset)

        last_error: Exception | None = None

        for url in dataset.urls:
            try:
                urlretrieve(url, cached_data_path)
                self._validate_checksum(dataset, cached_data_path)

                materialized = self._materialized(
                    dataset,
                    cached_data_path,
                )
                self.cache.write_manifest(dataset, materialized)
                return materialized

            except Exception as error:
                last_error = error

                if cached_data_path.exists():
                    cached_data_path.unlink()

        raise DatasetMaterializationError(
            f"Failed to download dataset {dataset.dataset_id}"
        ) from last_error


class SyntheticDatasetMaterializer(BaseDatasetMaterializer):
    def __init__(
        self,
        cache: DatasetCache,
        generators: SyntheticGeneratorRegistry,
    ) -> None:
        super().__init__(cache)
        self.generators = generators

    def materialize(
        self,
        dataset: BenchmarkDataset,
    ) -> MaterializedDataset:
        if not isinstance(dataset, SyntheticDataset):
            raise DatasetMaterializationError(
                "SyntheticDatasetMaterializer requires SyntheticDataset"
            )

        cached_data_path = self.cache.data_path_for(dataset)

        if self.cache.exists(dataset):
            self._validate_checksum(dataset, cached_data_path)
            return self._materialized(dataset, cached_data_path)

        generator = self.generators.get(dataset.generator)

        self.cache.prepare(dataset)

        generated_path = generator.generate(
            output_path=cached_data_path,
            parameters=dataset.parameters,
            seed=dataset.seed,
        )

        self._validate_checksum(dataset, generated_path)

        materialized = self._materialized(
            dataset,
            generated_path,
        )
        self.cache.write_manifest(dataset, materialized)

        return materialized


class DerivedDatasetMaterializer(BaseDatasetMaterializer):
    def __init__(
        self,
        cache: DatasetCache,
        dataset_registry: DatasetRegistry,
        materializer_registry: DatasetMaterializerRegistry,
        transformations: TransformationRegistry,
    ) -> None:
        super().__init__(cache)
        self.dataset_registry = dataset_registry
        self.materializer_registry = materializer_registry
        self.transformations = transformations

    def materialize(
        self,
        dataset: BenchmarkDataset,
    ) -> MaterializedDataset:
        if not isinstance(dataset, DerivedDataset):
            raise DatasetMaterializationError(
                "DerivedDatasetMaterializer requires DerivedDataset"
            )

        cached_data_path = self.cache.data_path_for(dataset)

        if self.cache.exists(dataset):
            self._validate_checksum(dataset, cached_data_path)
            return self._materialized(dataset, cached_data_path)

        parent = self.dataset_registry.get(dataset.parent_dataset_id)
        parent_materializer = self.materializer_registry.get(parent)
        parent_materialized = parent_materializer.materialize(parent)

        current_path = Path(parent_materialized.path)

        self.cache.prepare(dataset)

        for index, spec in enumerate(dataset.transformations):
            name = spec["name"]
            parameters = spec.get("parameters", {})

            transformation = self.transformations.get(name)

            is_last = index == len(dataset.transformations) - 1

            output_path = (
                cached_data_path
                if is_last
                else self.cache.path_for(dataset) / f"step_{index}"
            )

            current_path = transformation.apply(
                input_path=current_path,
                output_path=output_path,
                parameters=parameters,
            )

        if not dataset.transformations:
            copy_path(current_path, cached_data_path)
            current_path = cached_data_path

        self._validate_checksum(dataset, current_path)

        materialized = self._materialized(dataset, current_path)
        self.cache.write_manifest(dataset, materialized)

        return materialized


class DatasetMaterializerRegistry:
    def __init__(self) -> None:
        self._materializers: dict[type[BenchmarkDataset], DatasetMaterializer] = {}

    def register(
        self,
        dataset_class: type[BenchmarkDataset],
        materializer: DatasetMaterializer,
    ) -> None:
        self._materializers[dataset_class] = materializer

    def get(
        self,
        dataset: BenchmarkDataset,
    ) -> DatasetMaterializer:
        for dataset_class in type(dataset).__mro__:
            materializer = self._materializers.get(dataset_class)

            if materializer is not None:
                return materializer

        raise DatasetMaterializationError(
            f"No materializer registered for {type(dataset).__name__}"
        )
