from pathlib import Path

from maki.benchmark.data.models import Dataset
from maki.benchmark.workflow.base import BenchmarkPhase, BenchmarkWorkflow


class BuildGraphPhase(BenchmarkPhase):
    @property
    def name(self) -> str:
        return "build_graph"

    def execute(self, data: "Dataset"):
        forward, reverse = self._parse_files(data)
        self._cleanup(data)
            
    def _parse_files(self, data: "Dataset"):
        forward: Path | None = None
        reverse: Path | None = None
        
        for file in data.files:
            if file.role == "forward":
                forward = Path(file.url)
            elif file.role == "reverse":
                reverse = Path(file.url)
        
        assert forward, "Missing forward file"
        assert reverse, "Missing reverse file"
        
        print(forward, reverse)
        return forward, reverse
    
    def _cleanup(self, data: "Dataset"):
        if getattr(data, "cleanup_on_exit", False):
            cleaner = getattr(data, "cleanup", None)
            assert cleaner and callable(cleaner)
            cleaner()


class GraphBuildWorkflow(BenchmarkWorkflow):
    name = "graph_build"
    
    @staticmethod
    def compatible_data_types() -> list[str]:
        return ["paired-end"]

    @staticmethod
    def phases():
        return [
            BuildGraphPhase()
        ]
        