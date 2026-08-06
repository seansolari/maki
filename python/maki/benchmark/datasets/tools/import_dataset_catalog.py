
"""Allow dataset surveys to be maintained in a spreadsheet-like format and imported
directly into a versioned catalog manifest.
"""

from __future__ import annotations

import csv
from pathlib import Path

from maki.benchmark.datasets import (
    DatasetCatalog,
    DatasetMetadata,
    DatasetRegistry,
    DatasetType,
    RealDataset,
)


def load_table(path: Path):
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)

        for row in reader:
            yield row
            

def import_datasets(file: Path, output: Path, version: str):
    registry = DatasetRegistry()
  
    for row in load_table(file):
        dataset = RealDataset(
            dataset_id=row["dataset_id"],
            name=row["name"],
            version=row["version"],
            dataset_type=DatasetType.REAL,
            urls=[row["url"]],
            metadata=DatasetMetadata(
                description=row.get("description", ""),
                source=row.get("url"),
                checksum=row.get("checksum") or None,
                tags=[
                    tag.strip()
                    for tag in row.get("tags", "").split(";")
                    if tag.strip()
                ],
            ),
        )

        registry.register(dataset)

    catalog = DatasetCatalog(
        registry=registry,
        catalog_version=version,
    )

    catalog.save(output)
