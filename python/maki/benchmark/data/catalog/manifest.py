from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from ..models import Dataset, DatasetFile, Checksum


class ManifestError(Exception):
    pass


CORE_DATASET_KEYS = {
    "id",
    "name",
    "description",
    "domain",
    "data_type",
    "organism",
    "assay",
    "source",
    "license",
    "metadata",
    "files",
}


REQUIRED_MANIFEST_KEYS = ("schema_version", "datasets")
REQUIRED_DATASET_KEYS = ("id", "name", "domain", "data_type")
REQUIRED_FILE_KEYS = ("role", "url", "filename")
REQUIRED_CHECKSUM_KEYS = ("algorithm", "value")


def load_manifest(path: str | Path) -> list[Dataset]:
    """Load a JSON dataset manifest and return Dataset instances.

    The public API intentionally matches the previous YAML-backed loader:
    callers pass a path and receive ``list[Dataset]``. The input file is now
    expected to be JSON matching schema.json.
    """
    raw = _load_json(path)
    validate_manifest(raw)
    return [_parse_dataset(item) for item in raw["datasets"]]


def load_manifest_json(path: str | Path) -> list[Dataset]:
    """Explicit JSON-named alias for callers that prefer a descriptive name."""
    return load_manifest(path)


def _load_json(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    try:
        with path.open("r", encoding="utf-8") as f:
            raw = json.load(f)
    except json.JSONDecodeError as exc:
        raise ManifestError(f"Invalid JSON manifest: {exc.msg} at line {exc.lineno}, column {exc.colno}.") from exc
    except OSError as exc:
        raise ManifestError(f"Could not read manifest {path}: {exc}.") from exc

    if not isinstance(raw, dict):
        raise ManifestError("Manifest must be a JSON object.")
    return raw


def validate_manifest(raw: Mapping[str, Any]) -> None:
    """Validate the core manifest structure used by the parser.

    This lightweight validator enforces the same core contract captured in
    schema.json without requiring an additional runtime dependency. Projects that
    already use jsonschema can validate against schema.json before calling this
    loader for full JSON Schema error reporting.
    """
    if not isinstance(raw, Mapping):
        raise ManifestError("Manifest must be a JSON object.")

    for key in REQUIRED_MANIFEST_KEYS:
        if key not in raw:
            raise ManifestError(f"Manifest is missing {key}.")

    if not isinstance(raw["schema_version"], str) or not raw["schema_version"]:
        raise ManifestError("Manifest schema_version must be a non-empty string.")

    datasets = raw["datasets"]
    if not isinstance(datasets, list):
        raise ManifestError("Manifest datasets must be a list.")

    for index, item in enumerate(datasets):
        _validate_dataset(item, index)


def _validate_dataset(item: Any, index: int) -> None:
    if not isinstance(item, Mapping):
        raise ManifestError(f"Dataset at index {index} must be an object.")

    for key in REQUIRED_DATASET_KEYS:
        if key not in item:
            raise ManifestError(f"Dataset at index {index} is missing {key}.")
        if not isinstance(item[key], str) or not item[key]:
            raise ManifestError(f"Dataset at index {index} field {key} must be a non-empty string.")

    for optional_text_key in ("description", "organism", "assay", "source", "license"):
        if optional_text_key in item and item[optional_text_key] is not None and not isinstance(item[optional_text_key], str):
            raise ManifestError(f"Dataset {item.get('id', index)} field {optional_text_key} must be a string or null.")

    metadata = item.get("metadata", {})
    if metadata is not None and not isinstance(metadata, Mapping):
        raise ManifestError(f"Dataset {item.get('id', index)} metadata must be an object.")

    files = item.get("files", [])
    if not isinstance(files, list):
        raise ManifestError(f"Dataset {item.get('id', index)} files must be a list.")

    for file_index, file_item in enumerate(files):
        _validate_file(file_item, item.get("id", index), file_index)


def _validate_file(file_item: Any, dataset_id: str | int, file_index: int) -> None:
    if not isinstance(file_item, Mapping):
        raise ManifestError(f"File at index {file_index} for dataset {dataset_id} must be an object.")

    for key in REQUIRED_FILE_KEYS:
        if key not in file_item:
            raise ManifestError(f"File at index {file_index} for dataset {dataset_id} is missing {key}.")
        if not isinstance(file_item[key], str) or not file_item[key]:
            raise ManifestError(f"File field {key} at index {file_index} for dataset {dataset_id} must be a non-empty string.")

    size_bytes = file_item.get("size_bytes")
    if size_bytes is not None:
        if not isinstance(size_bytes, int) or size_bytes < 0:
            raise ManifestError(f"File size_bytes at index {file_index} for dataset {dataset_id} must be a non-negative integer or null.")

    checksum = file_item.get("checksum")
    if checksum is not None:
        if not isinstance(checksum, Mapping):
            raise ManifestError(f"Checksum at file index {file_index} for dataset {dataset_id} must be an object or null.")
        for key in REQUIRED_CHECKSUM_KEYS:
            if key not in checksum:
                raise ManifestError(f"Checksum at file index {file_index} for dataset {dataset_id} is missing {key}.")
            if not isinstance(checksum[key], str) or not checksum[key]:
                raise ManifestError(f"Checksum field {key} at file index {file_index} for dataset {dataset_id} must be a non-empty string.")


def _parse_dataset(item: Mapping[str, Any]) -> Dataset:
    files = [_parse_file(file_item) for file_item in item.get("files", [])]

    metadata = dict(item.get("metadata") or {})

    # Preserve dataset-specific top-level attributes without changing the
    # Dataset dataclass or the loader API. Consumers can read these from
    # dataset.metadata["attributes"].
    extra_attributes = {
        key: value for key, value in item.items() if key not in CORE_DATASET_KEYS
    }
    if extra_attributes:
        metadata.setdefault("attributes", {}).update(extra_attributes)

    return Dataset(
        id=item["id"],
        name=item["name"],
        description=item.get("description"),
        domain=item["domain"],
        data_type=item["data_type"],
        organism=item.get("organism"),
        assay=item.get("assay"),
        source=item.get("source"),
        license=item.get("license"),
        metadata=metadata,
        files=files,
    )


def _parse_file(file_item: Mapping[str, Any]) -> DatasetFile:
    checksum = None
    checksum_item = file_item.get("checksum")
    if checksum_item:
        checksum = Checksum(
            algorithm=checksum_item["algorithm"],
            value=checksum_item["value"],
        )

    return DatasetFile(
        role=file_item["role"],
        url=file_item["url"],
        filename=file_item["filename"],
        size_bytes=file_item.get("size_bytes"),
        checksum=checksum,
    )
