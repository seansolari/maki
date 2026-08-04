from maki.benchmark.workflow.base import BenchmarkPhase, BenchmarkWorkflow


class BuildGraphPhase(BenchmarkPhase):
    @property
    def name(self) -> str:
        return "build_graph"

    def execute(
        self,
        data,
    ):
        return (
            WeightedDeBruijnGraph
            .build(
                data.reads
            )
        )
        
        
class GraphBuildWorkflow(BenchmarkWorkflow):
    @property
    def name(self):
        return "graph_build"

    def phases(self):
        return [
            BuildGraphPhase()
        ]
        