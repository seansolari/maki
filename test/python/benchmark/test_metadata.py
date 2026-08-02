from maki.benchmark.metadata import collect_metadata


def test_metadata_populated():
    metadata = collect_metadata(threads=4)

    assert metadata.cpu_count > 0
    assert metadata.thread_count == 4
    assert metadata.python_version
