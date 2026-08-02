import json

from maki.benchmark import BenchmarkRunner
from maki.benchmark.stores.local import LocalResultStore


def test_result_file_written(tmp_path):
    runner = BenchmarkRunner()
    store = LocalResultStore(tmp_path)

    result = runner.run("example", lambda: None, store=store)

    output = tmp_path/ f"{result.run_id}.json"
    assert output.exists()

    data = json.loads(output.read_text())
    assert data["benchmark_name"] == "example"
