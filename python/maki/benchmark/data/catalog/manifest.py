# dataset_manager/manifest.py

from pathlib import Path
import yaml

from .models import Dataset, DatasetFile, Checksum


class ManifestError(Exception):
    pass


def load_manifest(path: str | Path) -> list[Dataset]:
    path = Path(path)

    with path.open("r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    if not isinstance(raw, dict):
        raise ManifestError("Manifest must be a mapping.")

    if "schema_version" not in raw:
        raise ManifestError("Manifest is missing schema_version.")

    if "datasets" not in raw:
        raise ManifestError("Manifest is missing datasets list.")

    datasets = []

    for item in raw["datasets"]:
        files = []

        for file_item in item.get("files", []):
            checksum = None
            if file_item.get("checksum"):
                checksum = Checksum(
                    algorithm=file_item["checksum"]["algorithm"],
                    value=file_item["checksum"]["value"],
                )

            files.append(
                DatasetFile(
                    role=file_item["role"],
                    url=file_item["url"],
                    filename=file_item["filename"],
                    size_bytes=file_item.get("size_bytes"),
                    checksum=checksum,
                )
            )

        datasets.append(
            Dataset(
                id=item["id"],
                name=item["name"],
                description=item.get("description"),
                domain=item["domain"],
                data_type=item["data_type"],
                organism=item.get("organism"),
                assay=item.get("assay"),
                source=item.get("source"),
                license=item.get("license"),
                metadata=item.get("metadata", {}),
                files=files,
            )
        )

    return datasets