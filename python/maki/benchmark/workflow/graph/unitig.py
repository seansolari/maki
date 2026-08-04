from maki.benchmark.workflow.base import BenchmarkPhase, BenchmarkWorkflow
from maki.benchmark.workflow.graph.build import BuildGraphPhase


class UnitigTransformPhase(BenchmarkPhase):
    @property
    def name(self):
        return "unitig_transform"

    def execute(self, data,):
        return UnitigView(data)
      

class UnitigWorkflow(BenchmarkWorkflow):
    @property
    def name(self):
        return "unitig_workflow"

    def phases(self):
        return [
            BuildGraphPhase(),
            UnitigTransformPhase(),
        ]
        