
from pathlib import Path
import urllib.request
import tempfile
import shutil

from .manifest import load_manifest
from .cache import DatasetCache
from .checksums import verify_checksum

"""
dm = DatasetManager(
    manifest_path="datasets.yaml",
    cache_dir="~/.cache/my_project/datasets",
    max_cache_size_gb=25,
)

for dataset in dm.list_datasets(data_type="shotgun_metagenomics"):
    print(dataset.id, dataset.name)

paths = dm.download("hmp_mock_shotgun_001")
print(paths)
"""


class DatasetNotFoundError(Exception):
    pass


class ChecksumMismatchError(Exception):
    pass


class DatasetManager:
    def __init__(
        self,
        manifest_path: str | Path,
        cache_dir: str | Path,
        max_cache_size_gb: float = 50,
    ):
        self.manifest_path = Path(manifest_path)
        self.datasets = load_manifest(self.manifest_path)
        self.by_id = {d.id: d for d in self.datasets}

        self.cache = DatasetCache(
            cache_dir=cache_dir,
            max_size_bytes=int(max_cache_size_gb * 1024**3),
        )

    def list_datasets(
        self,
        data_type: str | None = None,
        domain: str | None = None,
        organism: str | None = None,
        assay: str | None = None,
        tags: list[str] | None = None,
    ):
        results = self.datasets

        if data_type is not None:
            results = [d for d in results if d.data_type == data_type]

        if domain is not None:
            results = [d for d in results if d.domain == domain]

        if organism is not None:
            results = [d for d in results if d.organism == organism]

        if assay is not None:
            results = [d for d in results if d.assay == assay]

        if tags is not None:
            required = set(tags)
            results = [
                d for d in results
                if required.issubset(set(d.metadata.get("tags", [])))
            ]

        return results

    def get_dataset(self, dataset_id: str):
        try:
            return self.by_id[dataset_id]
        except KeyError:
            raise DatasetNotFoundError(f"Unknown dataset: {dataset_id}") from None

    def download(
        self,
        dataset_id: str,
        roles: list[str] | None = None,
        force: bool = False,
    ) -> dict[str, Path]:
        dataset = self.get_dataset(dataset_id)
        selected_files = dataset.files

        if roles is not None:
            role_set = set(roles)
            selected_files = [
                f for f in selected_files
                if f.role in role_set
            ]

        result = {}

        for file_record in selected_files:
            path = self.cache.path_for(dataset.id, file_record.filename)

            if path.exists() and not force:
                if self._is_valid_cached_file(path, file_record):
                    self.cache.touch(path)
                    result[file_record.role] = path
                    continue

                path.unlink()

            self._download_file(file_record.url, path)

            if not self._is_valid_cached_file(path, file_record):
                path.unlink(missing_ok=True)
                raise ChecksumMismatchError(
                    f"Checksum failed for {dataset.id}/{file_record.filename}"
                )

            self.cache.touch(path)
            result[file_record.role] = path

        self.cache.prune()
        return result

    def prune_cache(self) -> None:
        self.cache.prune()

    def clear_cache(self) -> None:
        self.cache.clear()

    def _is_valid_cached_file(self, path, file_record) -> bool:
        if file_record.size_bytes is not None:
            if path.stat().st_size != file_record.size_bytes:
                return False

        if file_record.checksum is not None:
            return verify_checksum(
                path,
                file_record.checksum.algorithm,
                file_record.checksum.value,
            )

        return True

    def _download_file(self, url: str, destination: Path) -> None:
        destination.parent.mkdir(parents=True, exist_ok=True)

        with tempfile.NamedTemporaryFile(
            dir=destination.parent,
            delete=False,
        ) as tmp:
            tmp_path = Path(tmp.name)

        try:
            urllib.request.urlretrieve(url, tmp_path)
            shutil.move(str(tmp_path), str(destination))
        finally:
            tmp_path.unlink(missing_ok=True)
