from typing import Optional

from maki.benchmark.workflow.base import BenchmarkPhase, BenchmarkWorkflow

from .build import BuildGraphPhase, validate_graph_build_extra_args


class UnitigTransformPhase(BenchmarkPhase):
    name = "unitig_transform"

    def transform(self, data):
        ...
      

class UnitigWorkflow(BenchmarkWorkflow):
    name = "unitig_workflow"
    
    phase_names = [
        BuildGraphPhase.name,
        UnitigTransformPhase.name
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
            BuildGraphPhase(self.k, self.s, self.threads, self.outdir),
            UnitigTransformPhase()
        ]

        