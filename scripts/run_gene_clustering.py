#!/usr/bin/env python3

"""Gene Clustering

Setup environment

```{console}
$ mamba create -n gene-clustering bioconda::mmseqs2 conda-forge::biopython
```
"""

from __future__ import annotations
from argparse import ArgumentParser
from collections import defaultdict
import enum
import gzip
from io import StringIO
from multiprocessing import Process, Queue
import multiprocessing
import os
import re
import queue
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tempfile import TemporaryDirectory
from typing import Dict, List, Mapping, Tuple


TRANSL_TABLE = 11

def file_exists(f: str) -> Tuple[str, bool]:
    return (f, os.path.exists(f))

class Gff3Splitter:
    def __init__(self, gff3_file: str, output_base: str):
        self.gff3_file = gff3_file

        tmp_prefix = os.path.join(output_base, os.path.basename(gff3_file))

        if tmp_prefix.endswith(".gz"):
            tmp_prefix = tmp_prefix[:-3]
            self.opener = lambda f: gzip.open(f, 'rt')
        else:
            self.opener = lambda f: open(f, 'rt')
        
        if tmp_prefix.endswith(".gff"):
            tmp_prefix = tmp_prefix[:-4]
        else:
            raise ValueError("Expected gff file, but %s is not .gff" % gff3_file)
        
        self.tmp_gff = tmp_prefix + "-tmp.gff"
        self.tmp_fna = tmp_prefix + "-tmp.fna"

    def __enter__(self):
        print("[Gff3Splitter] splitting %s into %s and %s" % (self.gff3_file, self.tmp_gff, self.tmp_fna))
        with self.opener(self.gff3_file) as ifh:
            with open(self.tmp_gff, 'wt') as ofh:
                for line in ifh:
                    if line != "##FASTA\n":
                        if not line.startswith("#"):
                            ofh.write(line)
                    else:
                        break
            with open(self.tmp_fna, 'wt') as ofh:
                for line in ifh:
                    ofh.write(line)
        
        return self.tmp_gff, self.tmp_fna

    def __exit__(self, *args, **kwargs):
        print("[Gff3Splitter] removing %s and %s" % (self.tmp_gff, self.tmp_fna))
        os.remove(self.tmp_gff)
        os.remove(self.tmp_fna)

def extract_protein_sequences_agat(gff_file: str, tempbase: str, logdir: str):
    with TemporaryDirectory(dir=tempbase) as workdir, Gff3Splitter(gff_file, workdir) as (input_gff_annots, input_gff_fna):
        prot_file = os.path.join(workdir, "proteins.fasta")
        subprocess.check_call(
            [
                "agat_sp_extract_sequences.pl",
                "-g", input_gff_annots,
                "-f", input_gff_fna,
                "-o", prot_file,
                "--table", "11",
                "-t", "cds",
                "-p"
            ],
            cwd=logdir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        protein_data = []
        with open(prot_file, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    seqid = re.findall(r"seq_id=(\S+)", line)[0]
                    line = f">{seqid}-{line[1:]}"
                protein_data.append(line)
        return protein_data

class AnnotList:
    def __init__(self) -> None:
        self.cds_by_parent = defaultdict(list)  # Parent -> list of parts
        self.meta_by_parent = {}                # Parent -> meta (strand, seqid, product, transl_table, locus_tag)

    @staticmethod
    def parse_attrs(attr_str):
        attrs = {}
        for kv in attr_str.strip().split(";"):
            if not kv:
                continue
            if "=" in kv:
                k, v = kv.split("=", 1)
                attrs[k] = v
        return attrs

    def add_annot(self, seqid: str, start: int, end: int, strand: str, phase: str, attrs: str):
        # Strand as +1/-1
        sgn = 1 if strand == "+" else -1 if strand == "-" else 0
        # Phase can be ".", "0", "1", "2"
        phase_val = None
        if phase in ("0", "1", "2"):
            phase_val = int(phase)

        a = self.parse_attrs(attrs)
        parent = a.get("Parent") or a.get("ID") or f"{seqid}:{start}-{end}"
        transl_table = int(a.get("transl_table", TRANSL_TABLE))
        product = a.get("product", "hypothetical protein")
        locus_tag = a.get("locus_tag") or a.get("protein_id") or parent

        # Store part
        self.cds_by_parent[parent].append({
            "seqid": seqid,
            "start": start,
            "end": end,
            "strand": sgn,
            "phase": phase_val,
            "attrs": a,
        })
        # Keep meta (consistent per parent)
        self.meta_by_parent[parent] = {
            "strand": sgn,
            "seqid": seqid,
            "product": product,
            "transl_table": transl_table,
            "locus_tag": locus_tag,
        }

    def yield_proteins(self, contigs: Dict[str, Seq]):
        for parent, parts in self.cds_by_parent.items():
            m = self.meta_by_parent[parent]
            seqid = m["seqid"]
            strand = m["strand"]
            transl_table = m["transl_table"]

            # Sort parts in transcription order
            if strand == 1:
                parts_sorted = sorted(parts, key=lambda x: x["start"])
            elif strand == -1:
                parts_sorted = sorted(parts, key=lambda x: x["start"], reverse=True)
            else:
                parts_sorted = sorted(parts, key=lambda x: x["start"])  # fallback

            cds_nt_chunks: List[Seq] = []

            for p in parts_sorted:
                # Extract raw slice (GFF is 1-based inclusive)
                nt = contigs[p["seqid"]][p["start"]-1 : p["end"]]

                # Apply phase trimming at the fragment's 5' end
                pv = p["phase"]
                if pv is not None and pv in (0, 1, 2) and pv != 0:
                    if strand == 1:
                        # Trim pv bases from left (5' is start on + strand)
                        nt = nt[pv:]
                    elif strand == -1:
                        # Trim pv bases from right (5' is end on - strand)
                        nt = nt[:-pv] if pv <= len(nt) else Seq("")  # avoid negative slice

                cds_nt_chunks.append(nt)

            # Concatenate and apply strand
            cds_nt = sum(cds_nt_chunks, Seq(""))
            if strand == -1:
                cds_nt = cds_nt.reverse_complement()

            # Translate
            try:
                aa = cds_nt.translate(table=transl_table, cds=True)
            except Exception:
                # Permissive fallback if CDS length not multiple of 3, or no terminal stop
                aa = cds_nt.translate(table=transl_table, to_stop=False)

            # Build record
            prot_id = f"{seqid}-{m['locus_tag']}"
            desc = f"product={m['product']}; Parent={parent}; contig={seqid}; n_parts={len(parts_sorted)}; strand={'+' if strand==1 else '-'}; transl_table={transl_table}"
            yield SeqRecord(aa, id=prot_id, description=desc)

def parse_gff_items(gff_file: str) -> AnnotList:
    # Collect annotations by Parent
    annots = AnnotList()
    with open(gff_file) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype != "CDS":
                continue
            start, end = int(start), int(end)
            annots.add_annot(seqid, start, end, strand, phase, attrs)
    return annots

def extract_protein_sequences_manual(gff_file: str, tempbase: str):
    with TemporaryDirectory(dir=tempbase) as workdir, Gff3Splitter(gff_file, workdir) as (input_gff_annots, input_gff_fna):
        # parse fasta and gff
        contigs = {rec.id: rec.seq for rec in SeqIO.parse(input_gff_fna, "fasta")}
        recs = parse_gff_items(input_gff_annots)
        # flush sequences to buffer
        buf = StringIO()
        for protein in recs.yield_proteins(contigs):
            SeqIO.write(protein, buf, "fasta")
        return buf.getvalue()

class ProteinRetriever(Process):
    def __init__(self, input_q: Queue, output_q: Queue, tmp_dir: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.input_q = input_q
        self.output_q = output_q
        self.tmp_dir = tmp_dir
        self.log_dir = os.path.join(self.tmp_dir, "logs", "agat")
        os.makedirs(self.log_dir, exist_ok=True)

    def run(self):
        while 1:
            try:
                input_gff = self.input_q.get(block=True, timeout=0.1)
            except queue.Empty:
                break
            protein_sequences = extract_protein_sequences_manual(input_gff, self.tmp_dir)
            self.output_q.put(protein_sequences, block=True, timeout=None)
        self.output_q.put(None, block=True, timeout=None)

class MMSeqs2RunConfiguration:
    def __init__(self, manifest_file: str, output_base: str, database_dir: str, num_procs: int):
        self.output_base = output_base

        self.temp_dir = os.path.join(self.output_base, "_mmseqs2_temp")
        self.search_dir = os.path.join(self.output_base, "search")
        self.cluster_dir = os.path.join(self.output_base, "cluster")
        self.classify_dir = os.path.join(self.output_base, "classify")

        self.dset_name = os.path.basename(output_base.rstrip("/"))
        self.mmseq_qry_database_dir = os.path.join(self.temp_dir, self.dset_name)
        self.mmseq_qry_database = os.path.join(self.mmseq_qry_database_dir, self.dset_name)

        self.database_dir = database_dir
        self.num_procs = num_procs

        # input data
        with open(manifest_file, 'r') as f:
            self.gff_files = [line.strip() for line in f]

        # check file existence
        with multiprocessing.Pool() as pool:
            fstats = pool.map(file_exists, self.gff_files)

        missing_files = [f for f, ex in fstats if not ex]
        if missing_files:
            missing_file_str = ", ".join(missing_files[:10])
            if len(missing_files) > 10:
                missing_file_str += ",..."
            print("[error] input gff files could not be found (%d): %s" % (len(missing_files), missing_file_str))
            exit(1)
        print("[MMSeqs2RunConfiguration] found %d input gff files" % len(self.gff_files))

        self.create_dirs()

    def create_dirs(self):
        if not os.path.exists(self.output_base):
            os.mkdir(self.output_base)
        if not os.path.exists(self.temp_dir):
            os.mkdir(self.temp_dir)
        if not os.path.exists(self.search_dir):
            os.mkdir(self.search_dir)
        if not os.path.exists(self.cluster_dir):
            os.mkdir(self.cluster_dir)
        if not os.path.exists(self.database_dir):
            os.mkdir(self.database_dir)
        if not os.path.exists(self.classify_dir):
            os.mkdir(self.classify_dir)

class MMSeqs2Database:
    def __init__(self, database_path: str):
        self.database_path = database_path

    def get_name(self) -> str:
        return os.path.basename(self.database_path.rstrip("/"))

    def list_headers(self) -> List[str] | None:
        lookup_file = self.database_path + ".lookup"

        if not os.path.exists(lookup_file):
            print("[error] could not find lookup file (expected at %s)" % lookup_file)
        else:
            with open(lookup_file, 'r') as f:
                return [line.strip().split()[1] for line in f]
            
    def exists(self) -> bool:
        return len(os.listdir(os.path.dirname(self.database_path))) > 0

class UniRefDatabase(enum.Enum):
    UniRef100 = 0
    UniRef90 = 1
    UniRef50 = 2

class MMSeqs2UniRefDatabase:
    def __init__(self, conf: MMSeqs2RunConfiguration, database: UniRefDatabase):
        if database == UniRefDatabase.UniRef100:
            self.database_name = "UniRef100"
            self.min_identity = 1.0
        elif database == UniRefDatabase.UniRef90:
            self.database_name = "UniRef90"
            self.min_identity = 0.9
        elif database == UniRefDatabase.UniRef50:
            self.database_name = "UniRef50"
            self.min_identity = 0.5
        else:
            raise ValueError("Unknown database type %s" % type(database))
        
        # membership criteria
        self.min_coverage = 0.8

        # file configuration
        self.database_base = os.path.join(conf.database_dir, self.database_name)
        self.database = os.path.join(self.database_base, self.database_name)
        self.linear_database = self.database + ".linidx"
        self.temp_dir = conf.temp_dir

        # run configuration
        self.thread_count = conf.num_procs

    def exists(self) -> bool:
        return os.path.exists(self.database_base)

    def download(self):
        if not os.path.exists(self.database_base):
            os.mkdir(self.database_base)

        cmd = [
            "mmseqs",
            "databases",
            self.database_name,
            self.database,
            self.temp_dir
        ]
        print("[MMSeqs2UniRefDatabase] running download command: %s" % " ".join(cmd))
        
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError as e:
            print("[MMSeqs2UniRefDatabase] intercepted exception, removing database dir %s" % self.database_base)
            os.rmdir(self.database_base)
            raise e
        
    def is_linearised(self):
        return os.path.exists(self.linear_database)

    def createlinindex(self):
        cmd = [
            "mmseqs",
            "createlinindex",
            self.database,
            self.temp_dir,
            "--threads", str(self.thread_count)
        ]
        print("[MMSeqs2UniRefDatabase] running indexing command: %s" % subprocess.list2cmdline(cmd))
        subprocess.check_output(cmd)

def make_MMSeqs2Database(conf: MMSeqs2RunConfiguration) -> MMSeqs2Database:
    if os.path.exists(conf.mmseq_qry_database_dir):
        print("[make_MMSeqs2Database] found database at %s" % conf.mmseq_qry_database_dir)
        return MMSeqs2Database(conf.mmseq_qry_database)
    
    os.mkdir(conf.mmseq_qry_database_dir)

    try:
        # input gff files

        input_gff_queue = Queue()
        for gff_file in conf.gff_files:
            input_gff_queue.put_nowait(gff_file)

        output_proteins_queue = Queue()

        # protein sequence extraction

        pool = [
            ProteinRetriever(input_gff_queue, output_proteins_queue, conf.temp_dir)
            for _ in range(conf.num_procs)
        ]

        p = subprocess.Popen(
            [
                "mmseqs", "createdb", "stdin", conf.mmseq_qry_database,
                "--dbtype", "1"
            ],
            stdin=subprocess.PIPE
        )
        
        with p.stdin:
            # start collecting protein sequences and piping to `MMSeqs2`
            
            for proc in pool:
                proc.start()

            num_finished_procs = 0

            # keep piping output until there are no remaining GFF files

            while not input_gff_queue.empty():
                seqs: str | None = output_proteins_queue.get(block=True, timeout=None)
                if seqs is None:
                    num_finished_procs += 1
                else:
                    p.stdin.write(seqs.encode("utf-8"))

            # input queue has been depleted, wait for no more output

            while (num_finished_procs < conf.num_procs) or (not output_proteins_queue.empty()):
                seqs: str | None = output_proteins_queue.get(block=True, timeout=None)
                if seqs is None:
                    num_finished_procs += 1
                else:
                    p.stdin.write(seqs.encode("utf-8"))

        for proc in pool:
            proc.join()

        ret_code = p.wait()
        if ret_code != 0:
            raise Exception("MMSeq2 createdb return non-zero exit code!")
        else:
            return MMSeqs2Database(conf.mmseq_qry_database)
    except Exception as e:
        os.rmdir(conf.mmseq_qry_database_dir)
        raise e

class MMSeqs2Result:
    def __init__(self, query: str, target: str, evalue: float, pident: float, alnlen: int, qcov: float, tcov: float):
        self.query = query
        self.target = target
        self.evalue = evalue
        self.pident = pident
        self.alnlen = alnlen
        self.qcov = qcov
        self.tcov = tcov

    @classmethod
    def from_tsv_line(cls, line: str) -> MMSeqs2Result:
        data = line.split("\t")
        assert len(data) == 7
        return MMSeqs2Result(data[0], data[1], float(data[2]), float(data[3]), int(data[4]), float(data[5]), float(data[6]))
    
    @classmethod
    def null_record(cls, accn: str) -> MMSeqs2Result:
        return MMSeqs2Result(accn, "unclassified", -1.0, -1.0, -1, -1.0, -1.0)

    def as_tuple(self):
        return (self.query, self.target, self.evalue, self.pident, self.alnlen, self.qcov, self.tcov)

    def to_string(self) -> str:
        return "%s\t%s\t%.4E\t%.3f\t%d\t%.3f\t%.3f\n" % (self.query, self.target, self.evalue, self.pident, self.alnlen, self.qcov, self.tcov)

class MMSeqs2ClusterDB:
    def __init__(self, input_db: str, min_identity: float, output_base: str):
        self.input_db = input_db
        self.input_db_name = os.path.basename(input_db)
        self.input_db_name_fmt = self.input_db_name if "-cluster-" not in self.input_db_name else \
            self.input_db_name.split("-cluster-")[0]

        self.output_name = self.input_db_name_fmt + ("-cluster-%.2f" % min_identity)
        self.output_base = os.path.join(output_base, self.output_name)
        self.clusterdb_path = os.path.join(self.output_base, self.output_name)
        self.tsv_file = self.clusterdb_path + "-results.tsv"

        self.min_identity = min_identity
        self.min_coverage = 0.8

        if not os.path.exists(self.output_base):
            os.mkdir(self.output_base)

    def has_results(self) -> bool:
        return len(os.listdir(self.output_base)) > 0
    
    def to_csv(self) -> Mapping[str, str]:
        if not os.path.exists(self.tsv_file):
            cmd = [
                "mmseqs", "createtsv",
                self.input_db, self.input_db, self.clusterdb_path, self.tsv_file
            ]
            print("[MMSeqs2ClusterDB] writing cluster results to %s" % self.tsv_file)
            subprocess.check_output(cmd)

        res = {}
        with open(self.tsv_file, 'r') as f:
            for line in f:
                rep, mem = line.strip().split("\t")
                res[mem] = rep
        return res
    
    def export_seeds(self, output_base: str) -> MMSeqs2Database:
        res_db_base = os.path.join(output_base, self.output_name)
        if not os.path.exists(res_db_base):
            os.mkdir(res_db_base)
        
        res_db_name = os.path.join(res_db_base, self.output_name)
        res = MMSeqs2Database(res_db_name)

        if not res.exists():
            cmd = [
                "mmseqs", "result2repseq",
                self.input_db, self.clusterdb_path, res.database_path
            ]
            print("[MMSeqs2ClusterDB] writing representative sequence database to %s" % res.database_path)
            subprocess.check_output(cmd)

        return res

class MMSeqs2SeachResultDB:
    def __init__(self, query_db: str, ref_db: str, output_base: str, output_name: str):
        self.query_db = query_db
        self.ref_db = ref_db
        self.output_base = os.path.join(output_base, output_name)
        self.results_path = os.path.join(self.output_base, output_name)
        self.formatted_results = self.results_path + "-output.txt"

        if not os.path.exists(self.output_base):
            os.mkdir(self.output_base)

    def has_results(self) -> bool:
        return len(os.listdir(self.output_base)) > 0
    
    def format_results(self) -> None:
        if not os.path.exists(self.formatted_results):
            cmd = [
                "mmseqs", "convertalis",
                self.query_db, self.ref_db, self.results_path, self.formatted_results,
                "--format-output", "query,target,evalue,pident,alnlen,qcov,tcov"
            ]
            print("[MMSeqs2SeachResultDB] running formatting command: %s" % subprocess.list2cmdline(cmd))
            subprocess.check_output(cmd)
        else:
            print("[MMSeqs2SeachResultDB] found formatted output at %s" % self.formatted_results)

class MMSeqs2:
    def __init__(self, conf: MMSeqs2RunConfiguration):
        self.search_dir = conf.search_dir
        self.cluster_dir = conf.cluster_dir
        self.temp_dir = conf.temp_dir
        self.thread_count = conf.num_procs

    def run_mmseqs_linclust(self, db: MMSeqs2Database, min_identity: float) -> MMSeqs2ClusterDB:
        res = MMSeqs2ClusterDB(db.database_path, min_identity, self.cluster_dir)

        if not res.has_results():
            cmd = [
                "mmseqs", "linclust",
                db.database_path, res.clusterdb_path, self.temp_dir,
                "--min-seq-id", str(res.min_identity),
                "-c", str(res.min_coverage),
                "--cov-mode", "1",
                "--kmer-per-seq", "81",
                "--threads", str(self.thread_count)
            ]
            print("[MMSeqs2] running MMSeqs2: %s" % subprocess.list2cmdline(cmd))
            subprocess.check_output(cmd)
        else:
            print("[MMSeqs2] cluster database already exists at %s" % res.clusterdb_path)
        
        return res

    def run_mmseqs_search(self, qry: MMSeqs2Database, ref: MMSeqs2UniRefDatabase) -> MMSeqs2SeachResultDB:
        results = MMSeqs2SeachResultDB(qry.database_path, ref.database, self.search_dir, ref.database_name + "-search")
        
        if not results.has_results():
            cmd = [
                "mmseqs", "linsearch", qry.database_path, ref.database,
                results.results_path, self.temp_dir,
                "--min-seq-id", str(ref.min_identity),
                "-c", str(ref.min_coverage),
                "--cov-mode", "1",
                "--kmer-per-seq", "81",
                "--threads", str(self.thread_count)
            ]
            print("[MMSeqs2] running MMSeqs2: %s" % " ".join(cmd))
            subprocess.check_output(cmd)
        else:
            print("[MMSeqs2] output already exists at %s" % results.results_path)

        return results

def take_highest_identity(results_file: str) -> Mapping[str, MMSeqs2Result]:
    dest = {}

    with open(results_file, "r") as f:
        for line in f:
            rec = MMSeqs2Result.from_tsv_line(line.strip())
            if rec.query not in dest:
                dest[rec.query] = rec
            elif rec.pident > dest[rec.query].pident:
                dest[rec.query] = rec

    return dest

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("-i", "--data", dest="manifest", type=str, required=True, help="Manifest of GFF files")
    parser.add_argument("--databases", dest="databases", type=str, required=True, help="Path to MMSeqs2 UniRef(100|90|50) databases folder")
    parser.add_argument("--threads", dest="num_parr", type=int, default=5, help="Number of processes")
    parser.add_argument("-o", "--output", dest="output", type=str, required=True, help="Output path")
    parser.add_argument("--cluster", action="store_true", help="Run clustering")
    parser.add_argument("--classify", action="store_true", help="Run classification")
    return parser.parse_args()

def MAIN_hierarchical_clustering(conf: MMSeqs2RunConfiguration, engine: MMSeqs2, qry_db: MMSeqs2Database):
    # Cluster at 90% AA idx

    qry_cls90 = engine.run_mmseqs_linclust(qry_db, 0.9)
    id2cls90 = qry_cls90.to_csv()
    qry_cls90_db = qry_cls90.export_seeds(conf.temp_dir)

    # Cluster seeds at 50% AA idx

    qry_cls50 = engine.run_mmseqs_linclust(qry_cls90_db, 0.5)
    seed2cls50 = qry_cls50.to_csv()

    # Write cluster results to disk

    run_name = qry_db.get_name()
    cls90_file = os.path.join(conf.cluster_dir, run_name + "-clusters-90aa.txt.gz")
    cls50_file = os.path.join(conf.cluster_dir, run_name + "-clusters-50aa.txt.gz")

    with gzip.open(cls90_file, 'wb') as f90, gzip.open(cls50_file, 'wb') as f50:
        for accn, rep90 in id2cls90.items():
            f90.write(b"%s\t%s\n" % (rep90.encode("utf-8"), accn.encode("utf-8")))
            rep50 = seed2cls50[rep90]
            f50.write(b"%s\t%s\n" % (rep50.encode("utf-8"), accn.encode("utf-8")))

def MAIN_search_uniref(conf: MMSeqs2RunConfiguration, engine: MMSeqs2, qry_db: MMSeqs2Database, uniref_databases: List[MMSeqs2UniRefDatabase]):
    qry_headers = qry_db.list_headers()
    if qry_headers is None:
        return
    
    # perform sequence search

    for ref_db in uniref_databases:
        # run search with MMSeqs2

        searchdb = engine.run_mmseqs_search(qry_db, ref_db)
        searchdb.format_results()

        # take highest identity hit per query

        best_hits = take_highest_identity(searchdb.formatted_results)

        # write results to output

        cluster_file = os.path.join(conf.classify_dir, ref_db.database_name + "-results.txt")
        assigned = 0
        with open(cluster_file, 'w') as f:
            for accn in qry_headers:
                try:
                    rec = best_hits[accn]
                    assigned += 1
                except KeyError:
                    rec = MMSeqs2Result.null_record(accn)
                f.write(rec.to_string())

        print("[main] wrote %d/%d assignments from database %s to %s" % (assigned, len(qry_headers), ref_db.database_name, cluster_file))

if __name__ == "__main__":
    args = parse_args()
    conf = MMSeqs2RunConfiguration(args.manifest, args.output, args.databases, args.num_parr)

    # make query database

    qry_db = make_MMSeqs2Database(conf)
    engine = MMSeqs2(conf)

    ## Cluster
    ## --------

    if (args.cluster):
        print("[main] running clustering...")
        MAIN_hierarchical_clustering(conf, engine, qry_db)

    ## Classify
    ## --------

    if (args.classify):
        print("[main] running classification...")

        # download reference databases
        uniref_databases: List[MMSeqs2UniRefDatabase] = []

        for unirefdb in (UniRefDatabase.UniRef50, UniRefDatabase.UniRef90):
            db = MMSeqs2UniRefDatabase(conf, unirefdb)
            if not db.exists():
                db.download()
            else:
                print("[main] found database %s at %s" % (db.database_name, db.database))
            if not db.is_linearised():
                db.createlinindex()
            uniref_databases.append(db)

        MAIN_search_uniref(conf, engine, qry_db, uniref_databases)
