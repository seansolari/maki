from __future__ import annotations

import time
import psutil


class MetricsCollector:
    def __init__(self) -> None:
        self._process = psutil.Process()

        self._start_wall = 0.0
        self._start_cpu = 0.0

    def start(self) -> None:
        self._start_wall = time.perf_counter()

        cpu = self._process.cpu_times()
        self._start_cpu = cpu.user + cpu.system

    def stop(self) -> tuple[float, float, int]:
        elapsed = time.perf_counter() - self._start_wall

        cpu = self._process.cpu_times()
        cpu_elapsed = (
            cpu.user + cpu.system
        ) - self._start_cpu

        peak_rss = self._process.memory_info().rss

        return elapsed, cpu_elapsed, peak_rss