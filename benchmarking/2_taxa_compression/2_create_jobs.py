#!/usr/bin/env python3
from argparse import ArgumentParser
from itertools import product
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("--manifest", type=str)
parser.add_argument("--out", type=str)
parser.add_argument("--env-init", type=str, default=None)
parser.add_argument("--email", type=str, default=None)
args = parser.parse_args()


# -------------------------
# OUTPUT
# -------------------------

out = Path(args.out)

jobs = out / "jobs"
jobs.mkdir(parents=True, exist_ok=True)

results = out / "results"
results.mkdir(parents=True, exist_ok=True)

slurm_out = out / "slurm"
slurm_out.mkdir(parents=True, exist_ok=True)


# -------------------------
# LOAD TEMPLATE
# -------------------------

with open(Path(__file__).parent / "job_template.sh", "r") as f:
    template = f.read()

if args.email:
    email_region = f"#SBATCH --mail-type=END,BEGIN,FAIL\n#SBATCH --mail-user={args.email}"
else:
    email_region = ""
    
    
# -------------------------
# GENERATE JOBS
# -------------------------

kmer_sizes = [11, 21, 31, 41]
ranks = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]

for k, rank in product(kmer_sizes, ranks):
    job_id = f"maki_{k}_{rank}"
    job_results = results / job_id
    job_results.mkdir(parents=True, exist_ok=True)
    
    script_content = template.format(
        job_name=job_id,
        slurm_out=slurm_out,
        email_region=email_region,
        threads=128,
        time="10:00:00",
        mem="128G",
        module_load=args.env_init or "",
        outdir=str(job_results),
        manifest=args.manifest,
        k=k,
        s=6,
        rank=rank,
        taxonomy="gtdb",
        releases="r214"
    )
    
    script_path = jobs / f"{job_id}.sh"

    with script_path.open("w") as f:
        f.write(script_content)

    script_path.chmod(0o755)
    