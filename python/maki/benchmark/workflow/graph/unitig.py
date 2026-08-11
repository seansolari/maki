from maki.benchmark.workflow.base import BenchmarkPhase, BenchmarkWorkflow

from .build import BuildGraphPhase


class UnitigTransformPhase(BenchmarkPhase):
    @property
    def name(self):
        return "unitig_transform"

    def execute(self, data,):
        ...
      

class UnitigWorkflow(BenchmarkWorkflow):
    name = "unitig_workflow"
    
    @staticmethod
    def compatible_data_types() -> list[str]:
        return ["paired-end"]

    @staticmethod
    def phases():
        return [
            BuildGraphPhase(),
            UnitigTransformPhase(),
        ]
        