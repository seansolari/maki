#!/usr/bin/env python3

from __future__ import annotations
from argparse import ArgumentParser
import builtins
from enum import Enum
from functools import partial
import gzip
import hashlib
import inspect
import json
import os
from pathlib import Path
import random
import re
import shutil
import string
import subprocess
import sys
import tempfile
from typing import Any, Callable, Dict, Generator, Iterable, List, Mapping, Tuple, TypeVar

# Run: ./run_benchmarking.py \
#   -B /mnt/DATASTORE/ssol0002/Projects/pam2s/build/bin/maki-main
#   -m /mnt/DATASTORE/GENOMES/GTDB/r220-reps/manifest.csv \
#   -g /mnt/DATASTORE/GENOMES/GTDB/r220-reps/manifest-gff.txt \
#   -o /mnt/DATASTORE/ssol0002/TEMP

GFF_EXTS = [".gff", ".gff3"]

"""Utils
"""

# OS

def yield_lines(file: str):
    with open(file, 'r') as f:
        for line in f:
            yield line.strip()

T = TypeVar('T')

def yield_from_each(*gens: Iterable[T]) -> Generator[Tuple[T, ...], None, None]:
    if len(gens) == 0:
        yield ()
    elif len(gens) == 1:
        for val in gens[0]:
            yield (val, )
    else:
        for val in gens[0]:
            for other_vals in yield_from_each(*gens[1:]):
                yield (val, ) + other_vals

def calculate_folder_size(folder: str) -> int:
    root_directory = Path(folder)
    return sum(f.stat().st_size for f in root_directory.glob('**/*') if f.is_file())

def get_md5(raw: str) -> str:
    return hashlib.md5(raw.encode('utf-8')).hexdigest()

def multinomial_diversify(p: List[int], n: int) -> List[int]:
    assert(n <= sum(p))

    q = [0 for _ in p]
    sampled = 0

    while sampled < n:
        qx = [i for i, (pv, qv) in enumerate(zip(p, q)) if qv < pv]
        for qi in random.choices(qx, k=n-sampled):
            if q[qi] < p[qi]:
                q[qi] += 1
                sampled += 1
    
    return q

class TemporaryDatasetFile:
    def __init__(self, file_name: str):
        self.file_path = os.path.join(tempfile.gettempdir(), file_name)

    def remove(self):
        print("[TemporaryDatasetFile] removing temporary dataset file %s" % self.file_path)
        os.remove(self.file_path)

"""Benchmarking

    - A `Benchmark` is comprised of a (command, dataset) pair, which is run on a parameter collection.
    - Datasets can be generated 'by hand', or by the results of other benchmarks.
    - Some parameters create new datasets, and some don't.
    - A `command`, combined with a `parameter set`, creates a `dataset`.
        - The `command` is just a functional that is called on a `parameter set`,
            it can be e.g. a Python class, subprocess call, etc.

"""

class ParameterOptions:
    def to_json(self):
        return {}
    
    @classmethod
    def from_json(cls, obj):
        return cls()

# Flag option for a parameter
class ParameterOptionFlag(ParameterOptions):
    def __init__(self, name: str):
        super().__init__()
        self.flag_name = name

    def to_json(self):
        return {"name": self.flag_name}
    
    @classmethod
    def from_json(cls, obj):
        return cls(obj["name"])

# All possible values for a paramter
class ParameterOptionValueSet(ParameterOptions):
    def __init__(self, name: str, values: List[Any]):
        super().__init__()
        self.value_name = name
        self.values = values

    def to_json(self):
        return {
            "name": self.value_name,
            "values": self.values
        }
    
    @classmethod
    def from_json(cls, obj):
        return cls(obj["name"], obj["values"])

# Instance from a parameter option set
class _ParameterOption:
    def __init__(self, parameter_name: str):
        self.name = parameter_name

# Parameter flag
class Flag(_ParameterOption):
    def __init__(self, parameter_name: str, flag_name: str):
        super().__init__(parameter_name)
        self.flag = flag_name

    def __repr__(self) -> str:
        return "ParameterValue(%s, %s)" % (self.name, self.flag)

# Value a parameter
class Value(_ParameterOption):
    def __init__(self, parameter_name: str, name: str, *value: Any):
        super().__init__(parameter_name)
        self.value_name = name
        self.value = value

    def __repr__(self) -> str:
        return "ParameterValue(%s, %s)" % (self.name, " ".join(str(v) for v in self.value))

class Parameter:
    def __init__(self, parameter_name: str, *args: ParameterOptions):
        self.name = parameter_name
        self.options = args

    def to_json(self):
        return {
            "parameter_name": self.name,
            "args": list(self.options)
        }
    
    @classmethod
    def from_json(cls, obj):
        return cls(obj["parameter_name"], *obj["args"])

    def __iter__(self):
        for opt in self.options:
            if isinstance(opt, ParameterOptionFlag):
                yield Flag(self.name, opt.flag_name)
            elif isinstance(opt, ParameterOptionValueSet):
                for opt_val in opt.values:
                    yield Value(self.name, opt.value_name, opt_val)
            else:
                raise TypeError("Unknown parameter type %s" % str(type(opt)))

class ParameterSet:
    def __init__(self, datum: List[Parameter] = [], neutral: List[Parameter] = []):
        self.dataset_generating_params = datum
        self.neutral_params = neutral

    def to_json(self):
        return {
            "datum": self.dataset_generating_params,
            "neutral": self.neutral_params
        }

    @classmethod
    def from_json(cls, obj):
        return cls(obj["datum"], obj["neutral"])

    def print_configuration(self):
        print("Parameter set:\n\tDataset generating parameters: %s\n\tNeutral tuning parameters: %s" % (
            ", ".join(p.name for p in self.dataset_generating_params),
            ", ".join(p.name for p in self.neutral_params)
        ))

# Benchmark turns a dataset into a parameter
class Dataset:
    def __init__(self, dataset_name: str, data: Any):
        self.name = dataset_name
        self.data = data

    def export_dataset(self) -> Value:
        raise NotImplementedError()
    
    def remove(self) -> None:
        raise NotImplementedError()

# A `Functional` takes as input a set of `_ParameterOption`s, runs an analysis, and returns output.
# This output may include:
#   - Benchmark statistics,
#   - `Dataset`.
# In this sense, a `Functional` behaves like a `ParameterOptionValueSet` that generates values 
# rather than stores them.
# These values are `Dataset`s.
# The type of analysis depends on the subclass of `Functional`.
class Functional:
    def __init__(self, output_parameter_name: str | None = None):
        self.output_parameter_name = output_parameter_name

    def to_json(self):
        return {"output_parameter_name": self.output_parameter_name}
    
    @classmethod
    def from_json(cls, obj):
        return cls(output_parameter_name = obj["output_parameter_name"])

    def print_configuration(self):
        if self.output_parameter_name is not None:
            print("Functional produces output parameter: %s" % self.output_parameter_name)
        else:
            print("Functional does not produce output parameter")

    def __call__(self, key_params: Tuple[_ParameterOption, ...], other_params: Tuple[_ParameterOption, ...]) -> Tuple[Any, FunctionalStats | None, Dataset | None]:
        key_param_dict = {p.name: p for p in key_params}
        other_param_dict = {p.name: p for p in other_params}
        try:
            raw_output, parsed_output, dataset = self.run(**key_param_dict, **other_param_dict)
            return raw_output, parsed_output, dataset
        except subprocess.CalledProcessError as e:
            raw_output = "output: %s\nstederr: %s\n" % (e.output.decode("utf-8"), e.stderr.decode("utf-8"))
            return raw_output, None, None

    def run(self, *args, **kwargs) -> Tuple[Any, FunctionalStats | None, Dataset | None]:
        raise NotImplementedError()

class FunctionalStats:
    def __init__(self, time: float | None, max_rss: int | None, cpu: float | None, disk_usage: float | None):
        self.time = time
        self.max_rss = max_rss
        self.cpu = cpu
        self.disk_usage = disk_usage

    def to_json(self):
        return {
            "time": self.time,
            "max_rss": self.max_rss,
            "cpu": self.cpu,
            "disk_usage": self.disk_usage
        }
    
    @classmethod
    def from_json(cls, obj):
        return cls(obj["time"], obj["max_rss"], obj["cpu"], obj["disk_usage"])
    
class BenchmarkResultKey:
    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self._key_order = sorted(self.kwargs.keys())

    def to_json(self):
        return {"kwargs": self.kwargs}
    
    @classmethod
    def from_json(cls, obj):
        return cls(**obj["kwargs"])

    def __key(self):
        return tuple(self.kwargs[k] for k in self._key_order)

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, rhs: Any):
        if isinstance(rhs, BenchmarkResultKey):
            return self.__key() == rhs.__key()
        return NotImplemented

class BenchmarkResults:
    def __init__(self, raw: Dict[BenchmarkResultKey, List[Any]] | None = None, parsed: Dict[BenchmarkResultKey, List[FunctionalStats]] | None = None):
        self.raw = {} if raw is None else raw
        self.parsed = {} if parsed is None else parsed

    def print_full(self):
        print("BenchmarkResults:\nRaw\n---")
        for key, rawlist in self.raw.items():
            print("\tkey: %s\n\t\t%s\n" % (str(key.kwargs), "\n\t\t".join(rawlist)))
        print("Parsed\n------")
        for key, statlist in self.parsed.items():
            print("\tkey: %s\n\t\t%s\n" % (str(key.kwargs), "\n\t\t".join(str(v.to_json()) for v in statlist)))

    def print(self):
        minraw, maxraw = None, None
        for k, v in self.raw.items():
            if minraw is None or len(v) < minraw:
                minraw = len(v)
            if maxraw is None or len(v) > maxraw:
                maxraw = len(v)
        if minraw is None:
            minraw = 0
        if maxraw is None:
            maxraw = 0
        
        minparsed, maxparsed = None, None
        for k, v in self.parsed.items():
            if minparsed is None or len(v) < minparsed:
                minparsed = len(v)
            if maxparsed is None or len(v) > maxparsed:
                maxparsed = len(v)
        if minparsed is None:
            minparsed = 0
        if maxparsed is None:
            maxparsed = 0

        print("Benchmark Results:\n\tRaw: %d (%d, %d)\n\tParsed: %d (%d, %d)" % (len(self.raw), minraw, maxraw, len(self.parsed), minparsed, maxparsed))

    def to_json(self):
        return {
            "raw": [[k, v] for k, v in self.raw.items()],
            "parsed": [[k, v] for k, v in self.parsed.items()]
        }
    
    @classmethod
    def from_json(cls, obj):
        return cls(
            raw = {k: v for (k, v) in obj["raw"]},
            parsed = {k: v for (k, v) in obj["parsed"]},
        )

    def insert_raw(self, key, value):
        try:
            self.raw[key].append(value)
        except KeyError:
            self.raw[key] = [value]

    def insert_parsed(self, key, value):
        try:
            self.parsed[key].append(value)
        except KeyError:
            self.parsed[key] = [value]

class Benchmark:
    def __init__(self, name: str, fn: Functional, ps: ParameterSet, reps: int = 3, output: BenchmarkResults | None = None, depends: List[Benchmark] | None = None):
        self.name = name

        # generators
        self.fn = fn
        self.params = ps

        # store output
        self.reps = reps
        self.output = BenchmarkResults() if output is None else output

        # dependent benchmarks
        self.dependent_benchmarks = [] if depends is None else depends

    def yield_benchmarks(self):
        yield self
        for bmark in self.dependent_benchmarks:
            yield from bmark.yield_benchmarks()

    def to_json(self):
        return {
            "name": self.name,
            "fn": self.fn,
            "ps": self.params,
            "reps": self.reps,
            "output": self.output,
            "depends": self.dependent_benchmarks
        }
    
    @classmethod
    def from_json(cls, obj):
        return cls(obj["name"], obj["fn"], obj["ps"], reps = obj["reps"], output = obj["output"], depends = obj["depends"])

    def print_configuration(self):
        print("BEGIN BENCHMARK --------------------------")
        print("Benchmark: %s (n=%d)" % (self.name, self.reps))
        self.fn.print_configuration()
        self.params.print_configuration()
        self.output.print()
        print("END BENCHMARK ----------------------------\n")

        for bmark in self.dependent_benchmarks:
            print("Benchmark dependent on %s:" % self.name)
            bmark.print_configuration()

    def add_benchmark(self, bmark: Benchmark):
        self.dependent_benchmarks.append(bmark)

    # Convert input dataset into a parameter
    def run(self, dataset_value: Value | None = None):
        runner = self._run_neutral_params if dataset_value is None else partial(self._run_neutral_params, dataset_value)
        for gen_params in yield_from_each(*self.params.dataset_generating_params):
            runner(*gen_params)
    
    # Neutral parameters are those that can vary while producing the same dataset
    def _run_neutral_params(self, *data_params: _ParameterOption):
        dataset = None

        for neu_params in yield_from_each(*self.params.neutral_params):
            result_key = self.make_key(*data_params, *neu_params)

            num_results = len(self.output.parsed.get(result_key, []))
            if num_results < self.reps:
                print("[Benchmark] found %d / %d results, running remaining replicates" % (num_results, self.reps))
                for i in range(self.reps - num_results):
                    raw_output, parsed_output, dataset = self.fn(data_params, other_params=neu_params)
                    self.output.insert_raw(result_key, raw_output)
                    self.output.insert_parsed(result_key, parsed_output)
            else:
                print("[Benchmark] found %d / %d results, no more replicates required" % (num_results, self.reps))
        
        if dataset is not None:
            try:
                dataset_param = dataset.export_dataset()
                for dependent_bmark in self.dependent_benchmarks:
                    dependent_bmark.run(dataset_param)
            finally:
                dataset.remove()

    @staticmethod
    def make_key(*args: _ParameterOption) -> BenchmarkResultKey:
        kd = {}

        for arg in args:
            if isinstance(arg, Flag):
                kd[arg.name] = arg.flag
            elif isinstance(arg, Value):
                kd[arg.name] = arg.value[0]
            else:
                raise TypeError("Unsupported result key type: %s" % type(arg))

        return BenchmarkResultKey(**kd)
    
    def yield_results(self):
        sorted_keys = None
        stat_attrs = ["time", "max_rss", "cpu", "disk_usage"]

        for key, results in self.output.parsed.items():
            if sorted_keys is None:
                sorted_keys = key._key_order
                yield sorted_keys + stat_attrs  # header
            key_comp = [str(key.kwargs[k]) for k in sorted_keys]
            for rep in results:
                if rep is not None:
                    yield key_comp + [str(getattr(rep, a)) for a in stat_attrs]

class CommandLineFunctional(Functional):
    def __init__(self, args: List[str], output_parameter_name: str | None = None):
        super().__init__(output_parameter_name)
        self.base_cmd = ["/usr/bin/time", "-f", "'%C\\t%e\\t%M\\t%P'"] + args

    def to_json(self):
        return {
            "args": self.base_cmd[3:],
            "output_parameter_name": self.output_parameter_name
        }

    @classmethod
    def from_json(cls, obj):
        return cls(obj["args"], output_parameter_name = obj["output_parameter_name"])

    def print_configuration(self):
        if self.output_parameter_name is not None:
            print("CL-functional runs command (%s) to produce output param: %s" % (" ".join(self.base_cmd), self.output_parameter_name))
        else:
            print("CL-functional runs command (%s)" % " ".join(self.base_cmd))

    @staticmethod
    def run_command(command: List[str]) -> str:
        return subprocess.check_output(command, stderr=subprocess.STDOUT).decode("utf-8")

    @staticmethod
    def parse_time_result(raw_output: str) -> Tuple[float, int, float]:
        time_line = raw_output.strip().split("\n")[-1]
        time_data = time_line.split("\t")
        
        if len(time_data) != 4:
            raise ValueError("Malformed /usr/bin/time output: %s" % time_line)

        time = float(time_data[1])
        rss = int(time_data[2])
        cpu = float(time_data[3][:-2])

        return time, rss, cpu

"""Functionals
"""

class DiversityClass(Enum):
    LOW = 1
    HIGH = 2

    def to_json(self):
        return {"diversity": str(self)}
    
    @classmethod
    def from_json(cls, obj):
        return getattr(cls, obj["diversity"].split(".")[1])

class TaxonomicRank(Enum):
    Species = 1
    Genus = 2
    Family = 3
    Order = 4
    Class = 5
    Phylum = 6

    def to_json(self):
        return {"rank": str(self)}
    
    @classmethod
    def from_json(cls, obj):
        return getattr(cls, obj["rank"].split(".")[1])

RANK2PREFIX = {
    TaxonomicRank.Species: "s__",
    TaxonomicRank.Genus: "g__",
    TaxonomicRank.Family: "f__",
    TaxonomicRank.Order: "o__",
    TaxonomicRank.Class: "c__",
    TaxonomicRank.Phylum: "p__"
}

class Genome:
    def __init__(self, name: str, file: str, taxonomy: str, meta):
        self.name = name
        self.file = file
        self.taxonomy = taxonomy
        self.metadata = meta

    def to_json(self):
        return {
            "name": self.name,
            "file": self.file,
            "taxonomy": self.taxonomy,
            "meta": self.metadata
        }
    
    @classmethod
    def from_json(cls, obj):
        return cls(obj["name"], obj["file"], obj["taxonomy"], obj["meta"])

    @staticmethod
    def file2name(file: str) -> str:
        fname = os.path.basename(file)

        if fname.endswith(".gz"):
            fname = fname[:-3]

        for ext in GFF_EXTS:
            if fname.endswith(ext):
                fname = fname[:-len(ext)]
                break
        else:
            raise Exception("File must be gff3 (%s)" % file)

        return fname

class GenomeDataset(Dataset):
    def __init__(self, dataset_name, diversity: DiversityClass, size: int, data: List[Genome]):
        super().__init__(dataset_name, data)
        self.dset_name = "%d-%s.txt" % (size, "lowdiv" if diversity == DiversityClass.LOW else "highdiv")
        self.temp_dset = TemporaryDatasetFile(self.dset_name)

    def export_dataset(self) -> Value:
        with open(self.temp_dset.file_path, "w") as f:
            for genome in self.data:
                f.write(genome.file + '\n')

        print("[GenomeSampleGroup] made temporary dataset file %s" % self.temp_dset.file_path)
        return Value(self.name, self.name, self.temp_dset.file_path)
    
    def remove(self) -> None:
        self.temp_dset.remove()

class GenomeSetGenerator(Functional):
    def __init__(self,
                 genomes_file: str | None = None,
                 metadata_file: str | None = None,
                 rank: TaxonomicRank = TaxonomicRank.Species,
                 meta_cols: List[str] | None = None,
                 groups: Mapping[str, List[Genome]] | None = None):
        super().__init__(output_parameter_name = "query_genomes")

        if groups is None:
            if genomes_file is None:
                raise ValueError("Must supply genome manifest")
            elif metadata_file is None:
                raise ValueError("Must supply genome metadata to sample from")
            
            print("[main] reading genomes from %s" % genomes_file)
            genome2file = { Genome.file2name(f): f for f in yield_lines(genomes_file) }

            print("[main] reading genome metadata from %s" % metadata_file)
            genomes = GenomeSetGenerator.read_metadata_for_genomes(metadata_file, genome2file, meta_cols=meta_cols)
            print("[main] retrieved metadata for %d genomes" % len(genomes))

            self.rank = rank
            self.groups = GenomeSetGenerator.group_by_taxonomy(genomes, rank)
        else:
            self.rank = rank
            self.groups = groups

    def to_json(self):
        return {"rank": self.rank, "groups": self.groups}
    
    @classmethod
    def from_json(cls, obj):
        return cls(rank = obj["rank"], groups = obj["groups"])

    """Initialisation Methods
    """

    @staticmethod
    def read_metadata_for_genomes(file: str, genome2file: Mapping[str, str], meta_cols: List[str] | None = None) -> List[Genome]:
        meta = []
        with open(file, 'r') as f:
            header = f.readline().strip().split(",")
            acc_col = header.index("Accession")
            tax_col = header.index("GTDB Taxonomy")
            
            mget: List[Tuple[str, int]] = []
            if meta_cols is not None:
                for meta_col in meta_cols:
                    try:
                        tp = (meta_col, header.index(meta_col))
                        mget.append(tp)
                    except ValueError:
                        pass

            for line in f:
                data = line.strip().split(",")
                acc = data[acc_col]
                tax = data[tax_col]
                if acc in genome2file:
                    meta.append(Genome(acc, genome2file[acc], tax, {mc: data[mci] for mc, mci in mget}))

        return meta
    
    @staticmethod
    def group_by_taxonomy(genomes: List[Genome], taxa: TaxonomicRank) -> Mapping[str, List[Genome]]:
        groups = {}

        if taxa not in RANK2PREFIX:
            raise ValueError("Unrecognised taxonomic rank: %s" % taxa)
        tregex = RANK2PREFIX[taxa] + "([^;]+)"
        print("[GenomeSetGenerator] collecting taxonomic group with regex: %s" % tregex)

        num_queries, num_assigned = len(genomes), 0

        while genomes:
            genome = genomes.pop()
            gp = re.findall(tregex, genome.taxonomy)
            if gp:
                try:
                    groups[gp[0]].append(genome)
                except KeyError:
                    groups[gp[0]] = [genome]
                num_assigned += 1

        print("[GenomeSetGenerator] assigned %d/%d to a taxonomic group" % (num_assigned, num_queries))
        return groups

    """Dataset Generators
    """

    def sample_low_diversity(self, n: int, seed: int) -> Generator[Genome, None, None]:
        random.seed(seed)

        # create groups that will be sampled

        gids = [gp for gp, genomes in self.groups.items() if len(genomes) > 0]
        gids.sort(key=lambda gid: len(self.groups[gid]), reverse=True)
        large_gids = [gid for gid in gids if len(self.groups[gid]) >= n]

        if sum(len(self.groups[gid]) for gid in gids) < n:
            raise ValueError("Not enough genomes to create group of size %d" % n)

        # select genomes using minimum diversity strategy

        num_selected = 0
        num_groups = 0

        while num_selected < n:
            # try sample large group
            
            if large_gids:
                gid = random.choice(large_gids)
                large_gids.remove(gid)
            else:
                gid = gids[0] # random.choice(gids)
            n_samples = min(n - num_selected, len(self.groups[gid]))

            yield from random.sample(self.groups[gid], n_samples)
            num_selected += n_samples
            num_groups += 1

            gids.remove(gid)

        print("[sample_low_diversity] sampled from %d groups" % num_groups)

    def sample_high_diversity(self, n: int, seed: int) -> Generator[Genome, None, None]:
        random.seed(seed)

        # create groups that will be sampled

        gids = [gp for gp, genomes in self.groups.items() if len(genomes) > 0]
        gids.sort(key = lambda gid: len(self.groups[gid]), reverse=True)

        if sum(len(self.groups[gid]) for gid in gids) < n:
            raise ValueError("Not enough genomes to create group of size %d" % n)

        for gii, gnum in enumerate(multinomial_diversify([len(self.groups[gid]) for gid in gids], n)):
            yield from random.sample(self.groups[gids[gii]], gnum)

    """API
    """

    def run(self, sample_diversity: Value, sample_size: Value) -> Tuple[List[Genome], None, GenomeDataset]:
        div, = sample_diversity.value
        size, = sample_size.value

        if div == DiversityClass.HIGH:
            sampler = self.sample_high_diversity
            seed = size
        elif div == DiversityClass.LOW:
            sampler = self.sample_low_diversity
            seed = ~size
        else:
            raise ValueError("Unrecognised diversity class: %s" % str(div))
        sampled_genomes = list(sampler(size, seed))
        return sampled_genomes, None, GenomeDataset(self.output_parameter_name, div, size, sampled_genomes)

class makiGraphDataset(Dataset):
    def export_dataset(self) -> Value:
        return Value(self.name, self.name, self.data)
    
    def remove(self) -> None:
        if not os.path.exists(self.data):
            print("[makiGraphDataset] WARNING: tried to remove graph at %s, but it could not be found" % self.data)
        else:
            shutil.rmtree(self.data)
            print("[makiGraphDataset] removed maki index at %s" % self.data)

class makiIndex(CommandLineFunctional):
    def __init__(self, output_base: str, command: List[str]):
        super().__init__(args=command, output_parameter_name="maki_graph")
        self.output_base = output_base
        self.command = command
    
    def to_json(self):
        return {
            "output_base": self.output_base,
            "command": self.command
        }
    
    @classmethod
    def from_json(cls, obj):
        return cls(obj["output_base"], obj["command"])

    @staticmethod
    def make_graph_name(query_file: str, kmer_size: str) -> str:
        query_name = os.path.basename(query_file)[:-4] # remove .txt suffix
        return "%s-k%s" % (query_name, kmer_size)

    @staticmethod
    def make_log_file(query_file: str, kmer_size: str, thread_count: str, build_mode: List[str]) -> str:
        query_name = os.path.basename(query_file)[:-4] # remove .txt suffix
        return "index-%s-k%s-t%s%s.log" % (query_name, kmer_size, thread_count, "".join(build_mode))

    def run(self, query_genomes: Value, kmer_size: Value, thread_count: Value, build_mode: Flag | Value):
        assert self.output_parameter_name is not None

        # parse input parameters

        qf, = query_genomes.value
        k = str(kmer_size.value[0])
        t = str(thread_count.value[0])

        if isinstance(build_mode, Flag):
            if build_mode.flag == "bulk":
                b = ["-b"]
            else:
                raise ValueError("Unrecognised flag name=%s, value=%s" % (build_mode.name, build_mode.flag))
        elif isinstance(build_mode, Value):
            b = ["-s", *map(str, build_mode.value)]

        # make dependent parameters

        graph_name = makiIndex.make_graph_name(qf, k)
        log_file = makiIndex.make_log_file(qf, k, t, b)
        
        # run index command

        graph_path = os.path.join(self.output_base, graph_name)
        log_path = os.path.join(self.output_base, log_file)

        cmd = [
            *self.base_cmd,
            "-q", qf,
            "-k", k,
            *b,
            "-t", t,
            "-o", graph_path,
            "-l", log_path
        ]
        print("[makiIndex] running command: %s" % " ".join(cmd))
        raw_output = self.run_command(cmd)
        time, rss, cpu = self.parse_time_result(raw_output)
        dspace = calculate_folder_size(graph_path)
        return raw_output, FunctionalStats(time, rss, cpu, dspace), makiGraphDataset(self.output_parameter_name, graph_path)

class makiCommand(CommandLineFunctional):
    def __init__(self, output_base: str, command: List[str]):
        super().__init__(args=command)
        self.output_base = output_base
        self.command = command

    def to_json(self):
        return {
            "output_base": self.output_base,
            "command": self.command
        }
    
    @classmethod
    def from_json(cls, obj):
        return cls(obj["output_base"], obj["command"])

    @staticmethod
    def make_output_file_name(*args, **kwargs) -> str:
        raise NotImplementedError()

    def run(self, maki_graph: Value, dist_mode: Flag , thread_count: Value) -> Tuple[str, FunctionalStats, None]:
        # parse input parameters

        graph, = maki_graph.value
        t = str(thread_count.value[0])
        d = []

        if dist_mode.flag == "genome":
            d.append("-W")
        elif dist_mode.flag != "gene":
            raise ValueError("Unrecognised value for `dist_mode`: %s" % repr(dist_mode))
        
        # make output file name

        output_dist = self.make_output_file_name(graph, d)
        output_path = os.path.join(self.output_base, output_dist)

        cmd = [
            *self.base_cmd,
            "-i", graph,
            "-o", output_path,
            "-t", t,
            *d
        ]
        print("[maki] running command: %s" % " ".join(cmd))

        raw_output = self.run_command(cmd)
        time, rss, cpu = self.parse_time_result(raw_output)
        return raw_output, FunctionalStats(time, rss, cpu, None), None
    
class makiDist(makiCommand):
    @staticmethod
    def make_output_file_name(graph_path: str, dist_mode: List[str]) -> str:
        graph_name = os.path.basename(graph_path).rstrip("/")
        return ("%s.dist.txt.gz" % graph_name) if not dist_mode else ("%s.dist.genome.txt.gz" % graph_name)

class makiCount(makiCommand):
    @staticmethod
    def make_output_file_name(graph_path: str, dist_mode: List[str]) -> str:
        graph_name = os.path.basename(graph_path).rstrip("/")
        return ("%s.count.txt.gz" % graph_name) if not dist_mode else ("%s.count.genome.txt.gz" % graph_name)

class makiCompress(makiCommand):
    @staticmethod
    def make_output_file_name(graph_path: str, dist_mode: List[str]) -> str:
        graph_name = os.path.basename(graph_path).rstrip("/")
        return ("%s.compress.txt.gz" % graph_name) if not dist_mode else ("%s.compress.genome.txt.gz" % graph_name)

class makiJaccard(CommandLineFunctional):
    def __init__(self, output_base: str, newick_file: str, command: List[str]):
        super().__init__(args=command, output_parameter_name="maki_jaccard")
        self.output_base = output_base
        self.newick_file = newick_file
        self.command = command

    def to_json(self):
        return {
            "output_base": self.output_base,
            "newick_file": self.newick_file,
            "command": self.command
        }
    
    @classmethod
    def from_json(cls, obj):
        return cls(obj["output_base"], obj["newick_file"], obj["command"])
    
    def get_genome_ids(self, graph_path: str) -> str:
        ids_file = os.path.join(graph_path, "genome_ids.txt")
        graph_name = os.path.basename(graph_path).rstrip("/")
        return shutil.copyfile(ids_file, os.path.join(self.output_base, graph_name + ".genome_ids.txt"))

    @staticmethod
    def make_output_file_name(graph_path: str) -> str:
        graph_name = os.path.basename(graph_path).rstrip("/")
        return f"{graph_name}.jaccard.txt.gz"

    def run(self, query_genomes: Value, kmer_size: Value, thread_count: Value):
        with tempfile.TemporaryDirectory(dir=self.output_base) as index_base:
            # create graph
            Index = makiIndex(index_base, [*self.command, "index", "-F"])
            _, _, graph = Index.run(query_genomes, kmer_size, thread_count, Value("build_mode", "suffix_size", 2))
            graph_path, = graph.export_dataset().value
            ids_file = self.get_genome_ids(graph_path)

            # count pairwise distances
            makiDist(self.output_base, [*self.command, "dist"])\
                .run(graph.export_dataset(), Flag("dist_mode", "genome"), thread_count)
            dist_file = os.path.join(self.output_base, makiDist.make_output_file_name(graph_path, ["-W"]))
            
            # count unique k-mers
            makiCount(self.output_base, [*self.command, "count"])\
                .run(graph.export_dataset(), Flag("dist_mode", "genome"), thread_count)
            count_file = os.path.join(self.output_base, makiCount.make_output_file_name(graph_path, ["-W"]))
            
            # calculate Jaccard indices
            k = str(kmer_size.value[0])
            jaccard_file = os.path.join(self.output_base, self.make_output_file_name(graph_path))

            cmd = [
                *self.base_cmd,
                "jaccard",
                "-i", ids_file,
                "-t", self.newick_file,
                "-c", count_file,
                "-d", dist_file,
                "-k", k,
                "-o", jaccard_file
            ]
            print(f"[maki] running command: {subprocess.list2cmdline(cmd)}")
            self.run_command(cmd)

            # delete graph files
            graph.remove()

        return None, None, None

class mashDist(CommandLineFunctional):
    def __init__(self, output_base, temp_base, command: List[str], thread_count: int):
        super().__init__(args=command)
        self.output_base = output_base
        self.temp_base = temp_base
        self.command = command
        self.thread_count = thread_count

    def to_json(self):
        return {
            "output_base": self.output_base,
            "temp_base": self.temp_base,
            "command": self.command,
            "thread_count": self.thread_count
        }
    
    @classmethod
    def from_json(cls, obj):
        return cls(obj["output_base"], obj["temp_base"], obj["command"], obj["thread_count"])
    
    @staticmethod
    def make_output_prefix(query_file: str, kmer_size: str) -> str:
        query_name = os.path.basename(query_file)[:-4] # remove .txt suffix
        return "%s-k%s" % (query_name, kmer_size)
    
    @staticmethod
    def gff2fasta(gff: str, outdir: str) -> str:
        name = os.path.basename(gff)
        if gff.endswith(".gz"):
            opener = lambda f: gzip.open(f, 'rb')
            name = name[:-3]
        else:
            opener = lambda f: open(f, 'rb')
        for gffext in GFF_EXTS:
            if name.endswith(gffext):
                name = name[:-len(gffext)]
                break
        else:
            raise ValueError("Require GFF file (%s)" % gff)
        
        output_file = os.path.join(outdir, name + ".fna")
        with open(output_file, 'wb') as o, opener(gff) as i:
            for line in i:
                if line == b"##FASTA\n":
                    break
            else:
                raise Exception("Could not find FASTA portion of gff input: %s" % gff)
            for line in i:
                o.write(line)
        return output_file

    def run(self, query_genomes: Value, kmer_size: Value):
        # parse input parameters

        qf, = query_genomes.value
        k = str(kmer_size.value[0])
        pref = mashDist.make_output_prefix(qf, k)

        with tempfile.TemporaryDirectory(dir=self.temp_base) as td:
            # create fasta files

            mash_manifest = []
            for gff_file in yield_lines(qf):
                mash_manifest.append(mashDist.gff2fasta(gff_file, td))
            manifest_file = os.path.join(td, "manifest.txt")
            with open(manifest_file, 'w') as f:
                f.write("\n".join(mash_manifest))

            # create sketch

            sketch_file = os.path.join(self.output_base, pref + ".msh")
            sketch_cmd = [
                *self.base_cmd, "sketch",
                "-p", str(self.thread_count),
                "-k", k,
                "-o", sketch_file,
                "-l", manifest_file
            ]
            sketch_output = self.run_command(sketch_cmd)
            sketch_time, sketch_rss, sketch_cpu = self.parse_time_result(sketch_output)

        # run mash dist

        dist_file = os.path.join(self.output_base, pref + ".dist.txt")
        dist_cmd = [
            *self.base_cmd, "dist",
            "-p", str(self.thread_count),
            sketch_file, sketch_file
        ]
        print("[mashDist] running command: %s" % subprocess.list2cmdline(dist_cmd))

        with open(dist_file, "w") as f:
            p = subprocess.run(dist_cmd, stdout=f, stderr=subprocess.PIPE)
        dist_output = p.stderr.decode("utf-8")
        dist_time, dist_rss, dist_cpu = self.parse_time_result(dist_output)
        
        # combine results

        output = "%s\n%s" % (sketch_output, dist_output)
        time = sketch_time + dist_time
        rss = max(sketch_rss, dist_rss)
        cpu = ((sketch_time * sketch_cpu) + (dist_time * dist_cpu)) / (sketch_time + dist_time)

        return output, FunctionalStats(time, rss, cpu, None), None

"""IO --------------------------------------------------------------------
"""

class ExtendedEncoder(json.JSONEncoder):
    def default(self, o):
        name = type(o).__name__
        try:
            encoder = getattr(o, "to_json")
        except AttributeError:
            return super().default(o)
        else:
            encoded = encoder()
            encoded["__extended_json_type__"] = name
            return encoded

class ExtendedDecoder(json.JSONDecoder):
    def __init__(self, **kwargs):
        kwargs["object_hook"] = self.object_hook
        super().__init__(**kwargs)

    def object_hook(self, obj):
        try:
            name = obj["__extended_json_type__"]
            decoder = getattr(getattr(sys.modules[__name__], name), "from_json")
        except (KeyError, AttributeError):
            return obj
        else:
            return decoder(obj)
        
def save_to_json(obj: Any, file: str) -> str:
    bytes_data = json.dumps(obj, ensure_ascii=False, indent=4, cls=ExtendedEncoder)\
        .encode("utf-8")
    
    if not file.endswith(".gz"):
        file += ".gz"

    with gzip.open(file, "wb") as f:
        f.write(bytes_data)
    print("[save_to_json] wrote data to %s" % file)

    return file

def load_from_json(file: str) -> Any:
    with gzip.open(file, "rb") as f:
        data = f.read().decode("utf-8")

    obj = json.loads(data, cls=ExtendedDecoder)
    print("[load_from_json] loaded data from %s" % file)
    
    return obj

"""=== Unit tests ===
"""

def TEST_EncodeDecode(seed: int):
    print("[TEST_EncodeDecode] running test with seed %d" % seed)
    random.seed(seed)

    # key that is used to group results
    
    attrs = {
        "int_att": lambda: random.randint(0, 100),
        "str_att": lambda: "".join(random.choices(string.ascii_letters, k = random.randint(0, 5)))
    }

    # randomly generate FunctionalStats Results
    
    generate_raw = lambda: "".join(random.choices(string.ascii_letters, k = random.randint(0, 10)))
    generate_stat = lambda: FunctionalStats(random.uniform(0.0, 100.0), random.randint(0, 100), random.uniform(0.0, 1000.0), random.uniform(0.0, 1e6))

    # replicates

    num_keys = 3
    reps = 3

    # generate fake data

    res = BenchmarkResults()

    for _ in range(num_keys):
        key = BenchmarkResultKey(**{k: g() for k, g in attrs.items()})

        for __ in range(reps):
            raw = generate_raw()
            stat = generate_stat()

            res.insert_raw(key, raw)
            res.insert_parsed(key, stat)

    # convert to JSON

    resdata = res.to_json()
    tar = BenchmarkResults.from_json(resdata)

    raweq = res.raw == tar.raw
    pareq = res.parsed == tar.parsed

    print(f"[TEST_EncodeDecode] Results...Raw equal:\t{raweq}")
    print(f"[TEST_EncodeDecode] Results...Parsed equal:\t{pareq}")

    if not raweq or not pareq:
        res.print_full()
        tar.print_full()

"""=== Dist Benchmark ===
"""

class DistOutputFileConfiguration:
    def __init__(self, base: str):
        self._base = base

        self._maki_base = os.path.join(self._base, "maki")
        self.maki_index = os.path.join(self._maki_base, "index")
        self.maki_dist = os.path.join(self._maki_base, "dist")
        self.maki_count = os.path.join(self._maki_base, "count")

        self.temp_base = os.path.join(self._base, "temp")

        self.groups_file = os.path.join(self._base, "groups.txt")
        self.json_file = os.path.join(self._base, "results.json.gz")
        self.stats_dir = self._base

        self.init_dirs()

    def init_dirs(self) -> None:
        for path in (self._base,
                     self._maki_base, self.maki_index, self.maki_dist, self.maki_count,
                     self.temp_base):
            if not os.path.exists(path):
                print("[DistOutputFileConfiguration] making output dir %s" % path)
                os.mkdir(path)

def ConfigureDistBenchmark(maki_binary: str, genomes_file: str, metadatafile: str, conf: DistOutputFileConfiguration, reps: int) -> Benchmark:
    # Group Creation Parameters
    gset_param_set = ParameterSet(
        datum=[
            Parameter("sample_diversity", ParameterOptionValueSet("div_class", [DiversityClass.LOW, DiversityClass.HIGH])),
            Parameter("sample_size", ParameterOptionValueSet("size", [100, 250, 500, 1000]))
        ])
    gset_fn = GenomeSetGenerator(genomes_file, metadatafile, rank=TaxonomicRank.Family)
    gset_generator = Benchmark("GenomeSetGenerator", gset_fn, gset_param_set, reps=1)

    # Global Index Parameters
    kmer_sizes = Parameter("kmer_size", ParameterOptionValueSet("k", [15, 21, 31]))
    thread_count_values = [5, 15, 25]
    thread_counts = Parameter("thread_count", ParameterOptionValueSet("threads", thread_count_values))

    # maki Index Parameters - depends on GenomeSetGenerator
    exix_generator = Benchmark("makiIndex",
                               makiIndex(conf.maki_index, [maki_binary, "index", "-F"]),
                               ParameterSet(datum=[kmer_sizes],
                                            neutral=[
                                                Parameter("build_mode",
                                                          ParameterOptionFlag("bulk"),
                                                          ParameterOptionValueSet("suffix_size", [1, 2])),
                                                thread_counts]),
                               reps=reps)
    gset_generator.add_benchmark(exix_generator)

    # maki Dist Parameters - depends on maki Index
    exdx_generator = Benchmark("makiDist",
                               makiDist(conf.maki_dist, [maki_binary, "dist"]),
                               ParameterSet(datum=[Parameter("dist_mode",
                                                             ParameterOptionFlag("genome"),
                                                             ParameterOptionFlag("gene"))],
                                            neutral=[thread_counts]),
                               reps=reps)
    exix_generator.add_benchmark(exdx_generator)

    # maki Count Parameters - depends on maki Index
    excx_generator = Benchmark("makiCount",
                               makiCount(conf.maki_count, [maki_binary, "count"]),
                               ParameterSet(datum=[Parameter("dist_mode",
                                                             ParameterOptionFlag("genome"),
                                                             ParameterOptionFlag("gene"))],
                                            neutral=[thread_counts]),
                               reps=reps)
    exix_generator.add_benchmark(excx_generator)
    
    return gset_generator

def dist(maki: str, genomes: str, metadata: str, outdir: str, reps: int):
    conf = DistOutputFileConfiguration(outdir)
    if os.path.exists(conf.json_file):
        bmark = load_from_json(conf.json_file)
    else:
        bmark = ConfigureDistBenchmark(maki, genomes, metadata, conf, reps)

    print("\nBenchmarking Configuration")
    print("==========================\n")
    bmark.print_configuration()
    bmark.run()

    # print genome set stats

    with open(conf.groups_file, "w") as f:
        f.write("sample_diversity\tsample_size\tgenome\ttaxonomy\n")
        for gkey, gsamples in bmark.output.raw.items():
            divstr = "highdiv" if gkey.kwargs["sample_diversity"] == DiversityClass.HIGH else "lowdiv"
            sizestr = str(gkey.kwargs["sample_size"])
            for genome in gsamples[0]:
                f.write("%s\t%s\t%s\t%s\n" % (divstr, sizestr, genome.name, genome.taxonomy))
    print("[main] wrote genome sample groups to %s" % conf.groups_file)
    
    # write results to csv

    for bm in bmark.yield_benchmarks():
        result_file_name = os.path.join(conf.stats_dir, bm.name + ".csv")
        with open(result_file_name, "w") as f:
            for result in bm.yield_results():
                f.write("\t".join(result) + "\n")
        print("[main] wrote stats from benchmark %s to %s" % (bm.name, result_file_name))

    # save raw data

    save_to_json(bmark, conf.json_file)

"""=== Compress Benchmark ===
"""

class CompressOutputFileConfiguration:
    def __init__(self, base: str):
        self._base = base

        self._maki_base = os.path.join(self._base, "maki")
        self.maki_index = os.path.join(self._maki_base, "index")
        self.maki_compress1 = os.path.join(self._maki_base, "compress1")

        self.temp_base = os.path.join(self._base, "temp")

        self.groups_file = os.path.join(self._base, "groups.txt")
        self.json_file = os.path.join(self._base, "results.json.gz")
        self.stats_dir = self._base

        self.init_dirs()

    def init_dirs(self) -> None:
        for path in (self._base,
                     self._maki_base, self.maki_index, self.maki_compress1,
                     self.temp_base):
            if not os.path.exists(path):
                print("[CompressOutputFileConfiguration] making output dir %s" % path)
                os.mkdir(path)

def ConfigureCompressBenchmark(maki_binary: str, genomes_file: str, metadatafile: str, conf: CompressOutputFileConfiguration, reps: int) -> Benchmark:
    # Group Creation Parameters
    gset_param_set = ParameterSet(
        datum=[
            Parameter("sample_diversity", ParameterOptionValueSet("div_class", [DiversityClass.LOW, DiversityClass.HIGH])),
            Parameter("sample_size", ParameterOptionValueSet("size", [100, 500, 1000]))
        ])
    gset_fn = GenomeSetGenerator(genomes_file, metadatafile, rank=TaxonomicRank.Family)
    gset_generator = Benchmark("GenomeSetGenerator", gset_fn, gset_param_set, reps=1)

    # Global Index Parameters
    kmer_sizes = Parameter("kmer_size", ParameterOptionValueSet("k", [15, 31]))
    thread_count_values = [32]
    thread_counts = Parameter("thread_count", ParameterOptionValueSet("threads", thread_count_values))

    # maki Index Parameters - depends on GenomeSetGenerator
    exix_generator = Benchmark("makiIndex",
                               makiIndex(conf.maki_index, [maki_binary, "index", "-F"]),
                               ParameterSet(datum=[kmer_sizes],
                                            neutral=[
                                                Parameter("build_mode",
                                                          ParameterOptionValueSet("suffix_size", [2])),
                                                thread_counts]),
                               reps=reps)
    gset_generator.add_benchmark(exix_generator)

    # maki Dist Parameters - depends on maki Index
    exdx_generator = Benchmark("makiCompress",
                               makiCompress(conf.maki_compress1, [maki_binary, "compress1"]),
                               ParameterSet(datum=[Parameter("dist_mode",
                                                             ParameterOptionFlag("gene"))],
                                            neutral=[thread_counts]),
                               reps=reps)
    exix_generator.add_benchmark(exdx_generator)
    
    return gset_generator

def compress1(maki: str, genomes: str, metadata: str, outdir: str, reps: int):
    conf = CompressOutputFileConfiguration(outdir)
    if os.path.exists(conf.json_file):
        bmark = load_from_json(conf.json_file)
    else:
        bmark = ConfigureCompressBenchmark(maki, genomes, metadata, conf, reps)

    print("\nBenchmarking Configuration")
    print("==========================\n")
    bmark.print_configuration()
    bmark.run()

    # print genome set stats

    with open(conf.groups_file, "w") as f:
        f.write("sample_diversity\tsample_size\tgenome\ttaxonomy\n")
        for gkey, gsamples in bmark.output.raw.items():
            divstr = "highdiv" if gkey.kwargs["sample_diversity"] == DiversityClass.HIGH else "lowdiv"
            sizestr = str(gkey.kwargs["sample_size"])
            for genome in gsamples[0]:
                f.write("%s\t%s\t%s\t%s\n" % (divstr, sizestr, genome.name, genome.taxonomy))
    print("[main] wrote genome sample groups to %s" % conf.groups_file)
    
    # write results to csv

    for bm in bmark.yield_benchmarks():
        result_file_name = os.path.join(conf.stats_dir, bm.name + ".csv")
        with open(result_file_name, "w") as f:
            for result in bm.yield_results():
                f.write("\t".join(result) + "\n")
        print("[main] wrote stats from benchmark %s to %s" % (bm.name, result_file_name))

    # save raw data

    save_to_json(bmark, conf.json_file)

"""=== Mash Benchmark ===
"""

class MashOutputFileConfiguration:
    def __init__(self, base: str):
        self._base = base

        self.maki_base = os.path.join(self._base, "maki")
        self.mash_base = os.path.join(self._base, "mash")
        self.temp_base = os.path.join(self._base, "temp")

        self.groups_file = os.path.join(self._base, "groups.txt")
        self.json_file = os.path.join(self._base, "results.json.gz")
        self.stats_dir = self._base

        self.init_dirs()

    def init_dirs(self) -> None:
        for path in (self._base, self.maki_base, self.mash_base, self.temp_base):
            if not os.path.exists(path):
                print("[MashOutputFileConfiguration] making output dir %s" % path)
                os.mkdir(path)

def ConfigureMashBenchmark(maki_binary: str, mash_binary: str, genomes_file: str, metadatafile: str, newick: str, conf: MashOutputFileConfiguration) -> Benchmark:
    # Group Creation Parameters
    gset_param_set = ParameterSet(
        datum=[
            Parameter("sample_diversity", ParameterOptionValueSet("div_class", [DiversityClass.LOW, DiversityClass.HIGH])),
            Parameter("sample_size", ParameterOptionValueSet("size", [1000]))
        ])
    gset_fn = GenomeSetGenerator(genomes_file, metadatafile, rank=TaxonomicRank.Family)
    gset_generator = Benchmark("GenomeSetGenerator", gset_fn, gset_param_set, reps=1)

    # Global Index Parameters
    kmer_sizes = Parameter("kmer_size", ParameterOptionValueSet("k", [15, 21, 31]))
    thread_counts = Parameter("thread_count", ParameterOptionValueSet("threads", [35]))

    # maki distances
    maki_generate = Benchmark("makiJaccard",
                              makiJaccard(conf.maki_base, newick, [maki_binary]),
                              ParameterSet(datum=[kmer_sizes],
                                           neutral=[thread_counts]),
                              reps=1)
    gset_generator.add_benchmark(maki_generate)

    # mash distances
    mash_generator = Benchmark("mashDist",
                               mashDist(conf.mash_base, conf.temp_base, [mash_binary], 35),
                               ParameterSet(datum=[kmer_sizes]),
                               reps=1)
    gset_generator.add_benchmark(mash_generator)
    
    return gset_generator

def mash(maki: str, mash: str, genomes: str, metadata: str, newick: str, outdir: str):
    conf = MashOutputFileConfiguration(outdir)
    if os.path.exists(conf.json_file):
        bmark = load_from_json(conf.json_file)
    else:
        bmark = ConfigureMashBenchmark(maki, mash, genomes, metadata, newick, conf)

    print("\nBenchmarking Configuration")
    print("==========================\n")
    bmark.print_configuration()
    bmark.run()

    # print genome set stats

    with open(conf.groups_file, "w") as f:
        f.write("sample_diversity\tsample_size\tgenome\ttaxonomy\n")
        for gkey, gsamples in bmark.output.raw.items():
            divstr = "highdiv" if gkey.kwargs["sample_diversity"] == DiversityClass.HIGH else "lowdiv"
            sizestr = str(gkey.kwargs["sample_size"])
            for genome in gsamples[0]:
                f.write("%s\t%s\t%s\t%s\n" % (divstr, sizestr, genome.name, genome.taxonomy))
    print("[main] wrote genome sample groups to %s" % conf.groups_file)
    
    # write results to csv

    for bm in bmark.yield_benchmarks():
        result_file_name = os.path.join(conf.stats_dir, bm.name + ".csv")
        with open(result_file_name, "w") as f:
            for result in bm.yield_results():
                f.write("\t".join(result) + "\n")
        print("[main] wrote stats from benchmark %s to %s" % (bm.name, result_file_name))

    # save raw data

    save_to_json(bmark, conf.json_file)

""" Main ----------------------------------------------------------------------
"""

def run_CLI(*modules: Callable):
    # define parser
    parser = ArgumentParser()
    subcommand = parser.add_subparsers(dest="module", help="Benchmarking module to run", required=True)
    for module in modules:
        subparser = subcommand.add_parser(module.__name__.lower())
        for param in inspect.signature(module).parameters.values():
            is_required = param.kind == inspect.Parameter.POSITIONAL_OR_KEYWORD and param.default is inspect.Parameter.empty
            subparser.add_argument("--" + param.name,
                                   dest=param.name,
                                   type=getattr(builtins, param.annotation, type(None)),
                                   required=is_required)
    
    # call parser
    args = parser.parse_args()
    for module in modules:
        if args.module == module.__name__.lower():
            print(f"[main] Running module: {args.module}")
            fwd_args = {}
            for param in inspect.signature(module).parameters.values():
                val = getattr(args, param.name)
                if val is not None:
                    fwd_args[param.name] = val
            module(**fwd_args)
            break
    else:
        parser.print_help()
        exit(1)

if __name__ == "__main__":
    run_CLI(TEST_EncodeDecode,
            dist,
            compress1,
            mash)
