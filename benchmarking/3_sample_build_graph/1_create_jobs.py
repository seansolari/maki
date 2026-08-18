#!/usr/bin/env python3
from argparse import ArgumentParser
from itertools import product
import math
import os
from pathlib import Path
import subprocess
from typing import List


parser = ArgumentParser()
parser.add_argument("--datasets", nargs='+', type=str)
parser.add_argument("--out", type=str)
parser.add_argument("--env-init", type=str, default=None)
parser.add_argument("--limit", type=int, default=-1)
parser.add_argument("--email", type=str, default=None)
parser.add_argument("--no-run", action="store_true")
parser.add_argument("--force", action="store_true")
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
# HELPER FUNCTIONS
# -------------------------

def yield_sample_pairs(files: List[str]):
    files.sort()
    assert len(files)%2 == 0
    for i in range(len(files)//2):
        fwd = files[(2*i)]
        rev = files[(2*i)+1]
        sample_name = os.path.basename(fwd).split(".")[0]
        yield sample_name, fwd, rev

def pair_samples(files: List[str]):
    return list(yield_sample_pairs(files))


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

samples = pair_samples(args.datasets)
kmer_sizes = [21, 31, 41]
threads = [32, 64]
suffix_sizes = [4, 5, 6, 7, 8]
read_counts: List[int | None] = [100_000, 500_000, 1_000_000, None]

job_count = 0
submitted_count = 0
submitted_limit = args.limit if args.limit > -1 else math.inf

for i, ((sname, fwd, rev), k, thread, s, rcount) in enumerate(product(samples, kmer_sizes, threads, suffix_sizes, read_counts)):
    if submitted_count >= submitted_limit:
        break
    
    job_id = f"maki-bmark-{i}"
    job_results = results / job_id

    # Skip completed runs
    
    if (job_results / "done.flag").exists():
        stat = (job_results / "done.flag").read_text().strip()
        if (stat == "submitted"):
            if args.force:
                print(f"Re-submitting job {job_id}.")
            else:
                print(f"Skipping {job_id} (already submitted, to override add `--force`)")
                continue
        else:
            assert stat == "complete"
            print(f"Skipping {job_id} (already done)")
            continue

    job_results.mkdir(parents=True, exist_ok=True)

    script_content = template.format(
        job_name=job_id,
        slurm_out=slurm_out,
        email_region=email_region,
        threads=thread,
        outdir=str(job_results),
        module_load=args.env_init or "",
        forward=fwd,
        reverse=rev,
        subsample_region="" if rcount is None else f"--subsample {rcount}",
        k=k,
        s=s,
        subsample=rcount or "All"
    )

    script_path = jobs / f"{job_id}.sh"
    
    with script_path.open("w") as f:
        f.write(script_content)

    script_path.chmod(0o755)

    # Submit job
    if not args.no_run:
        subprocess.run(["sbatch", script_path])
        (job_results / "done.flag").write_text("submitted")

        submitted_count += 1
    
    job_count += 1

print(f"Generated jobs: {job_count}")
print(f"Submitted jobs: {submitted_count}")
