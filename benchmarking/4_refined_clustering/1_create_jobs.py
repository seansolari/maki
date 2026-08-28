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

kmer_sizes = [21, 31, 41]
ani_thresholds = [0.95, 0.96, 0.97, 0.98, 0.99]

for k, ani in product(kmer_sizes, ani_thresholds):
    job_id = f"maki_{k}_{ani}"
    job_results = results / job_id
    job_results.mkdir(parents=True, exist_ok=True)
    
    script_content = template.format(
        job_name=job_id,
        slurm_out=slurm_out,
        email_region=email_region,
        threads=128,
        time="24:00:00",
        mem="128G",
        module_load=args.env_init or "",
        outdir=str(job_results),
        manifest=args.manifest,
        k=k,
        taxonomy="gtdb",
        releases="r214",
        ani=ani,
        s=6,
    )
    
    script_path = jobs / f"{job_id}.sh"

    with script_path.open("w") as f:
        f.write(script_content)

    script_path.chmod(0o755)
    