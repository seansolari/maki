from maki.benchmark.workflow.base import BenchmarkPhase, BenchmarkWorkflow


class BuildGraphPhase(BenchmarkPhase):
    @property
    def name(self) -> str:
        return "build_graph"

    def execute(
        self,
        data,
    ):
        ...
        
        
class GraphBuildWorkflow(BenchmarkWorkflow):
    @property
    def compatible_data_types(self) -> list[str]:
        return ["shotgun_metagenomic"]
    
    @property
    def name(self):
        return "graph_build"

    def phases(self):
        return [
            BuildGraphPhase()
        ]
        