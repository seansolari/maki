from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Iterator, Literal, overload


# ---------------------------------------------------------------------------
# Filename parser interface
# ---------------------------------------------------------------------------

@dataclass(slots=True)
class FastqFileHandle:
    sample_name: str
    lane: str | None
    fragment: str | None
    suffix: str
    path: Path


@dataclass(slots=True)
class OrientedFastqFileHandle(FastqFileHandle):
    orientation: int


class FastqAttributeParsingException(Exception):
    pass


class FastqOrientationException(Exception):
    pass


@overload
def parse_fastq_attributes(file: Path, force_oriented: Literal[False]) -> FastqFileHandle: ...

@overload
def parse_fastq_attributes(file: Path, force_oriented: Literal[True]) -> OrientedFastqFileHandle: ...

def parse_fastq_attributes(file: Path, force_oriented: bool = False) -> FastqFileHandle | OrientedFastqFileHandle:
    """
    Supplied by the application.

    It should parse a sequence filename and return a dataclass-like object
    exposing:

        sample_name
        lane
        orientation
        fragment
        suffix

    Replace this stub with the actual implementation or import.
    """
    matches = re.match(
        r"^(?P<sample_name>.+?)(?:[_.-](?P<lane>L\d\d\d))?[_.-](?P<orientation>R?[12])?"
        r"(?:[_.-](?P<fragment>\d\d\d))?(?P<suffix>\.(?:f(?:ast)?q)(?:\.gz)?)$",
        file.name
    )
    
    if not matches:
        raise FastqAttributeParsingException(f"Could not parse data from fastq name: {file.name}")
    
    vals = matches.groupdict()
    
    if force_oriented and ("orientation" not in vals):
        raise FastqOrientationException(f"Could not parse sequence orientation (R1, R2) from fastq name: {file.name}")
    
    if "orientation" in vals:
        return OrientedFastqFileHandle(
            vals["sample_name"],
            vals.get("lane", None),
            vals.get("fragment", None),
            vals["suffix"],
            file,
            int(vals["orientation"][-1]),
        )
    else:
        return FastqFileHandle(
            vals["sample_name"],
            vals.get("lane", None),
            vals.get("fragment", None),
            vals["suffix"],
            file
        )
    
    
# ----------------------------------------------------------------------
# Data groups
# ----------------------------------------------------------------------

@dataclass(slots=True)
class SequenceFile:
    suffix: str
    path: Path


@dataclass(slots=True)
class PairedReads:
    r1: SequenceFile
    r2: SequenceFile


@dataclass(slots=True)
class LaneFragment:
    fragment: str
    paired_data: PairedReads


@dataclass(slots=True)
class ExperimentLane:
    lane: str
    fragments: dict[str, LaneFragment] = field(default_factory=dict)


@dataclass(slots=True)
class SampleData:
    sample_name: str
    lanes: dict[str, ExperimentLane] = field(default_factory=dict)
    
    def iter(self) -> Iterator[tuple[str, str, PairedReads]]:
        for lane in self.lanes.values():
            for fragment in lane.fragments.values():
                yield lane.lane, fragment.fragment, fragment.paired_data
                        
    def flatten(self) -> Iterator[PairedReads]:
        for lane in self.lanes.values():
            for fragment in lane.fragments.values():
                yield fragment.paired_data


# ----------------------------------------------------------------------
# Parsing
# ----------------------------------------------------------------------

def discover_paired_fastqs(directory: Path) -> list[OrientedFastqFileHandle]:
    files: list[OrientedFastqFileHandle] = []
    
    for file in directory.rglob("*"):
        try:
            files.append(parse_fastq_attributes(file, force_oriented=True))
        except FastqAttributeParsingException:
            pass
        
    return files


def parse_paired_fastq_list(manifest_path: Path) -> list[OrientedFastqFileHandle]:
    files: list[OrientedFastqFileHandle] = []
    
    with manifest_path.open("rt") as f:
        for line in f:
            file = Path(line.strip())
            
            if not file.exists():
                raise OSError(f"Input file {file} not found.")
            
            files.append(parse_fastq_attributes(file, force_oriented=True))
        
    return files


class MultipleUnspecifiedFragmentsException(Exception):
    pass


class InvalidOrientationsException(Exception):
    pass


def group_fastqs(handles: list[OrientedFastqFileHandle]) -> list[SampleData]:
    # Group by fragment, lane and sample
    
    grouped: dict[
        tuple[str, str, str],
        dict[int, list[OrientedFastqFileHandle]]
    ] = defaultdict(dict)
    
    for handle in handles:
        key = (
            handle.sample_name,
            handle.lane or "null",
            handle.fragment or "null"
        )
        grouped[key][handle.orientation].append(handle)
    
    # Ensure paired
    
    sample_data: dict[str, SampleData] = {}
    
    for (sample, lane, fragment), fragment_data in grouped.items():
        orientations: list[int] = []
        
        for orientation, files in fragment_data.items():
            if len(files) > 1:
                raise MultipleUnspecifiedFragmentsException(
                    f"Multiple fragments ({len(files)}) for sample={sample}, "
                    f"lane={lane}, fragment={fragment}, orientation={orientation}."
                )
                
            orientations.append(orientation)
        
        orientations.sort()
        if orientations != [1, 2]:
            raise InvalidOrientationsException(
                f"Invalid orientations ({'/'.join(map(str, orientations))}, "
                f"expect 1/2) for sample={sample}, lane={lane}, "
                f"fragment={fragment}."
            )
        
        sample_data.setdefault(sample, SampleData(sample))\
            .lanes.setdefault(lane, ExperimentLane(lane))\
            .fragments[fragment] = LaneFragment(
                fragment,
                PairedReads(
                    r1=SequenceFile(
                        fragment_data[1][0].suffix,
                        fragment_data[1][0].path
                    ),
                    r2=SequenceFile(
                        fragment_data[2][0].suffix,
                        fragment_data[2][0].path
                    )
                )
            )

    return list(sample_data.values())
