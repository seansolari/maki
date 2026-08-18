
from dataclasses import dataclass
from pathlib import Path
import re


@dataclass(slots=True)
class SampleNameComponents:
    sample_name: str
    lane: str | None
    orientation: int
    fragment: str | None
    suffix: str


def fastq_sample_name(file: Path):
    matches = re.match(
        r"^(?P<sample_name>.+?)(?:[_.-](?P<lane>L\d\d\d))?[_.-](?P<orientation>R?[12])"
        r"(?:[_.-](?P<fragment>\d\d\d))?(?P<suffix>\.(?:f(?:ast)?q)(?:\.gz)?)$",
        file.name
    )
    
    if not matches:
        raise RuntimeError(f"Could not parse data from fastq name: {file.name}")
    
    vals = matches.groupdict()
    
    return SampleNameComponents(
        vals["sample_name"],
        vals.get("lane", None),
        int(vals["orientation"][-1]),
        vals.get("fragment", None),
        vals["suffix"]
    )
