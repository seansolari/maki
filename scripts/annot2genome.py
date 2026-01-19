#!/usr/bin/env python3

from argparse import ArgumentParser
import gzip
import multiprocessing
import os
import re

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

def yield_gzip(file: str):
    with gzip.open(file, 'rb') as f:
        for line in f:
            yield line.decode("utf-8").strip()

def yield_txt(file: str):
    with open(file, 'r') as f:
        for line in f:
            yield line.strip()

def annt2genome(file: str):
    # collect data

    genome_name = file2genomename(file)

    print("[annt2genome] processing %s (%s)" % (genome_name, file))

    annot2accn = {}
    accn2genome = {}

    if file.endswith(".gz"):
        line_gen = yield_gzip(file)
    else:
        line_gen = yield_txt(file)

    for line in line_gen:
        if line == "##FASTA":
            break
        gff = line.split("\t")
        accn = gff[0]
        annotid = re.findall(r"ID=([^;\n]+)", gff[-1])[0]
        annot2accn[annotid] = accn
        accn2genome[accn] = genome_name

    # encode data

    data = [
        ("%s\t%s\t%s\n" % (annotid, accn, accn2genome[accn])).encode("utf-8")
        for annotid, accn in annot2accn.items()
    ]

    return data

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("-f", "--manifest", dest="manifest", type=str, required=True, help="File listing input genomes")
    parser.add_argument("-o", "--output", dest="output", type=str, required=True, help="Where to save output")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # input genomes

    with open(args.manifest, 'r') as f:
        genomes = [line.strip() for line in f]

    # configure output folder

    output_file: str = args.output

    if not output_file.endswith(".gz"):
        output_file += '.gz'

    with multiprocessing.Pool() as pool, gzip.open(output_file, "wb") as f:
        for dset in pool.imap(annt2genome, genomes):
            for line in dset:
                f.write(line)
    
    print("Wrote data to %s" % output_file)

if __name__ == "__main__":
    main()
