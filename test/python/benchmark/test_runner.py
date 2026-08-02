from maki.benchmark import BenchmarkRunner


def test_run_without_persistence():
    runner = BenchmarkRunner()

    result = runner.run("example", lambda: sum(range(1000)))

    assert result.benchmark_name == "example"
    assert result.metrics.wall_time_seconds >= 0
