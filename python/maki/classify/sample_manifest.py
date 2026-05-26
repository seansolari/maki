import csv
from dataclasses import dataclass


@dataclass(frozen=True)
class Sample:
    name: str
    forward: str
    reverse: str


class SampleManifest:
    def __init__(self, samples):
        self.samples = samples

    @classmethod
    def from_csv(cls, path):
        samples = []
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                samples.append(
                    Sample(
                        name=row["sample"],
                        forward=row["forward"],
                        reverse=row["reverse"]
                    )
                )
        return cls(samples)
