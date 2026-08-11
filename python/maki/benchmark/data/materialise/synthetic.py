
from dataclasses import dataclass

from maki.benchmark.data.models import DatasetFile, DatasetReference
from maki.benchmark.data.synthetic.dataset import SyntheticBenchmarkDataset
from maki.benchmark.data.synthetic.sequencing_registry import SequencingSimulatorOptions, build_simulator
from maki.benchmark.data.synthetic.topologies import SyntheticSequence

from .materialiser import DatasetMaterialiser


@dataclass
class SyntheticDatasetToken(DatasetReference):
    sequences: list[SyntheticSequence]
    sequencing_opts: SequencingSimulatorOptions


class SyntheticDatasetMaterialiser(DatasetMaterialiser):
    input_type = SyntheticDatasetToken
    
    def __init__(self, cleanup_on_exit: bool) -> None:
        super().__init__()
        self.cleanup_on_exit = cleanup_on_exit
    
    def materialise(self, token: SyntheticDatasetToken) -> SyntheticBenchmarkDataset:
        simulator = build_simulator(token.sequencing_opts)
        dataset = simulator.generate(token.sequences)
        
        return SyntheticBenchmarkDataset(
            token.id,
            token.name,
            token.description,
            token.domain,
            dataset.data_type,
            None,
            None,
            None,
            None,
            metadata={
                "topologies": sorted({s.topology for s in token.sequences}),
                "sequence_count": len(token.sequences),
            },
            files=[DatasetFile(role, str(file), file.name) for role, file in dataset.files.items()],
            sequences=token.sequences,
            cleanup_on_exit=self.cleanup_on_exit,
            temporary_directory=dataset.temporary_directory,
        )
