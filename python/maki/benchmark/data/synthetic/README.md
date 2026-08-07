# Synthetic benchmark generation

This lightweight submodule generates synthetic linear sequences that stress canonical graph structures and passes them through a minimal sequencing layer.

## Example

```python
from maki.benchmark.data.synthetic import (
    TopologyParameters,
    ShortPairedEndSimulator,
    generate_synthetic_benchmark,
)

params = TopologyParameters(k=31, sequence_count=2, length=2000, seed=42)
simulator = ShortPairedEndSimulator(read_length=150, read_count=1000, insert_size=350)

with generate_synthetic_benchmark(
    simulator=simulator,
    topology_parameters=params,
    cleanup_on_exit=True,
) as dataset:
    assert dataset.data_type == "shotgun_metagenome"
    workflow_input_files = dataset.files
    # {'forward': Path(...), 'reverse': Path(...)}
```

Supported topology generators:

- `linear_chain`
- `bubble`
- `deep_branching`
- `high_degree_repeat`
- `random`

Supported sequencing simulators:

- `ShortPairedEndSimulator`, returns `{'forward': ..., 'reverse': ...}`
- `ShortUnpairedSimulator`, returns `{'unpaired': ...}`
- `LongReadSimulator`, returns `{'long_reads': ...}`
