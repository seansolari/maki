from maki.benchmark.workflow.base import BenchmarkPhase, BenchmarkWorkflow
from maki.benchmark.workflow.graph.build import BuildGraphPhase


class UnitigTransformPhase(BenchmarkPhase):
    @property
    def name(self):
        return "unitig_transform"

    def execute(self, data,):
        ...
      

class UnitigWorkflow(BenchmarkWorkflow):
    @property
    def compatible_data_types(self) -> list[str]:
        return ["shotgun_metagenomic"]
    
    @property
    def name(self):
        return "unitig_workflow"

    def phases(self):
        return [
            BuildGraphPhase(),
            UnitigTransformPhase(),
        ]
        