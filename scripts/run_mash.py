#!/usr/bin/env python3

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from functools import partial
import gzip
from multiprocessing import Pool
import os
import shutil
import subprocess
from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import List, Mapping

GFF_EXTS = [".gff", ".gff3"]

def file2genomename(file: str) -> str:
    fname = os.path.basename(file)

    if fname.endswith(".gz"):
        fname = fname[:-3]

    for ext in GFF_EXTS:
        if fname.endswith(ext):
            fname = fname[:-len(ext)]
            break
    else:
        raise ValueError("Could not detect gff file from name: %s" % file)
    
    return fname

def mash2name(file: str) -> str:
    return os.path.basename(os.path.splitext(file)[0])

def check_mash() -> bool:
    return shutil.which("mash") is not None

def read_manifest(manifest: str) -> List[str]:
    with open(manifest, 'r') as f:
        return [line.strip() for line in f]
    
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

def mash_sketch(sample_name: str, query: str, k: int, s: int, threads: int, outdir: str) -> str:
    prefix = os.path.join(outdir, "%s.k%d.s%d" % (sample_name, k, s))
    ofile = prefix + ".msh"
    if not os.path.exists(ofile):
        cmd = [
            "mash", "sketch",
            "-p", str(threads),
            "-l",
            "-o", prefix,
            "-k", str(k),
            "-s", str(s),
            query
        ]
        subprocess.check_call(cmd)
    return ofile

def mash_dist(sketch_file: str, threads: int, outdir: str) -> str:
    sample_name = os.path.basename(sketch_file)[:-4]
    gzofile = os.path.join(outdir, sample_name + ".txt.gz")

    if not os.path.exists(gzofile):
        # run mash dist
        cmd = [
            "mash", "dist",
            "-p", str(threads),
            sketch_file, sketch_file
        ]

        bytes_output = subprocess.check_output(cmd)

        with gzip.open(gzofile, 'wb') as f:
            f.write(bytes_output)
    
    return gzofile

def mash_to_table(dist_file: str, genome2id: Mapping[str, int]) -> str:
    import numpy as np

    # output location

    table_file = dist_file[:-7] + ".tbl.txt"

    # prepare table data

    num_genomes = len(genome2id)
    table = np.ones((num_genomes, num_genomes), dtype=np.float64)

    # populate table 

    assert dist_file.endswith(".gz")
    with gzip.open(dist_file, 'rb') as f:
        for line in f:
            qry, tgt, dst, pval, hshstr = line.decode("utf-8").strip().split("\t")
            r, c = genome2id[mash2name(qry)], genome2id[mash2name(tgt)]
            if c < r:
                r, c = c, r
            table[r, c] = float(dst)

    # write table to disk

    id2genome = {i: v for v, i in genome2id.items()}

    with open(table_file, "w") as f:
        f.write("%d\n" % num_genomes)
        for i in range(num_genomes):
            f.write("%s\t%s\n" % (id2genome[i], "\t".join(str(v) for v in table[i])))

    return table_file

def check_rapidnj() -> bool:
    return shutil.which("rapidnj") is not None

def _format_rapidnj_out(nwk: str) -> str:
    seq_types = [".fna", ".faa"]
    comp_types = [".tar.gz", ".gz", ""]

    nwk = nwk.replace("'", "")

    for seq_type in seq_types:
        for comp_type in comp_types:
            ext = seq_type + comp_type + ":"
            nwk = nwk.replace(ext, ":")

    return nwk

def rapidnj(dist_file: str, threads: int, outdir: str) -> str:
    sample_name = os.path.basename(dist_file)[:-8]
    ofile = os.path.join(outdir, sample_name + ".nwk")
    if not os.path.exists(ofile):
        cmd = [
            "rapidnj",
            dist_file,
            "-i", "pd",
            "-o", "t",
            "-c", str(threads)
        ]
        nwk = _format_rapidnj_out(subprocess.check_output(cmd).decode("utf-8"))

        with open(ofile, 'w') as f:
            f.write(nwk)

    return ofile

def gzip_file(file: str) -> None:
    subprocess.check_call(["gzip", file])

def parse_args():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--manifest", dest="manifest", type=str, required=True, help="File listing input genomes")
    parser.add_argument("-k", dest="k", type=int, default=21, help="K-mer size")
    parser.add_argument("-s", dest="s", type=int, default=10000, help="Sketch size")
    parser.add_argument("-o", "--out", dest="out", type=str, required=True, help="Output directory")
    parser.add_argument("-t", "--threads", dest="threads", type=int, required=True, help="Number of threads")
    parser.add_argument("--tree", action="store_true", dest="tree", default=False, help="Run RapidNJ")
    return parser.parse_args()

def main():
    if not check_mash():
        print("[main] ERROR: Could not find Mash installation")
        exit(1)

    args = parse_args()

    if not os.path.exists(args.out):
        os.mkdir(args.out)

    # load genomes

    genomes = read_manifest(args.manifest)

    # write fasta files to temp

    with TemporaryDirectory(prefix=os.getcwd()) as tmpdir:
        with Pool(processes=args.threads) as pool:
            queries = list(pool.imap_unordered(partial(gff2fasta, outdir=tmpdir), genomes))

        # write list of input fasta files

        with NamedTemporaryFile(mode="w+", suffix=".txt") as tmpfile:
            for query in queries:
                tmpfile.write(query + "\n")
            tmpfile.flush()

            # create mash sketch

            sample_name: str = os.path.basename(args.manifest)
            if sample_name.endswith(".txt"):
                sample_name = sample_name[:-4]
            sketch_file = mash_sketch(sample_name, tmpfile.name, args.k, args.s, args.threads, args.out)
            print("[main] created sketch file: %s" % sketch_file)

    # run mash dist

    dist_file = mash_dist(sketch_file, args.threads, args.out)
    print("[main] created distance table: %s" % dist_file)

    if args.tree:
        genome2id = {x: i for i, x in enumerate(sorted([file2genomename(f) for f in genomes]))}

        # convert mash format to table
        table_file = mash_to_table(dist_file, genome2id)

        # run rapidnj
        if not check_rapidnj():
            print("[main] ERROR: Could not find RapidNJ installation")
            exit(1)

        nwk_file = rapidnj(table_file, args.threads, args.out)
        subprocess.check_call(["gzip", table_file])
        print("[main] created tree: %s" % nwk_file)

if __name__ == "__main__":
    main()
