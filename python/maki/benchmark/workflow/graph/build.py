import os
from typing import Optional

from maki.benchmark.data.models import Dataset, DatasetFile
from maki.benchmark.workflow.base import BenchmarkPhase, BenchmarkWorkflow


def validate_graph_build_extra_args(args: dict):
    parsed = {}
    missing_args = []
    
    for arg, arg_type in (("k", int), ("s", int)):
        try:
            parsed[arg] = arg_type(args[arg])
        except KeyError:
            missing_args.append((arg, arg_type.__name__))
            
    for arg, arg_type in (("threads", int), ("graph_out", str)):
        try:
            parsed[arg] = arg_type(args[arg])
        except KeyError:
            pass
    
    extra_args = [arg for arg in args if arg not in parsed]
    
    return parsed, missing_args, extra_args


class BuildGraphPhase(BenchmarkPhase):
    name = "build_graph"
    
    def __init__(self, k: int, s: int, threads: int, outdir: Optional[str]) -> None:
        super().__init__()
        
        self.k = k
        self.s = s
        self.threads = threads
        self.outdir = outdir

    def transform(self, data: "Dataset"):
        forward, reverse = self._parse_files(data.files)
        outdir = self.outdir or os.path.dirname(forward)
        graph_output = os.path.join(outdir, data.id.replace(' ', ''))
        
        cmd = [
            "maki", "sample-to-debruijn",
            forward, reverse,
            str(self.k), str(self.s), str(self.threads), graph_output
        ]
        
        return cmd, graph_output
    
    def _parse_files(self, files: list[DatasetFile]):
        forward: str | None = None
        reverse: str | None = None
        
        for file in files:
            if file.role == "forward":
                forward = str(file.url)
            elif file.role == "reverse":
                reverse = str(file.url)
        
        assert forward, "Missing forward file"
        assert reverse, "Missing reverse file"
        
        return forward, reverse


class GraphBuildWorkflow(BenchmarkWorkflow):
    name = "graph_build"
    
    phase_names = [
        BuildGraphPhase.name
    ]
    
    compatible_data_types = [
        "paired-end"
    ]
    
    @staticmethod
    def validate_args(args: dict):
        return validate_graph_build_extra_args(args)
    
    def __init__(self, k: int, s: int, threads: int = 1, graph_out: Optional[str] = None) -> None:
        super().__init__()
        
        self.k = k
        self.s = s
        self.threads = threads
        self.outdir = graph_out

    def phases(self):
        return [
            BuildGraphPhase(self.k, self.s, self.threads, self.outdir)
        ]
        