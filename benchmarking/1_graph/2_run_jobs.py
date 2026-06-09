#!/usr/bin/env python3
from argparse import ArgumentParser
from dataclasses import dataclass
import os
from pathlib import Path
import subprocess
from typing import List


parser = ArgumentParser()
parser.add_argument("--config", type=str)
parser.add_argument("--out", type=str)
parser.add_argument("--env-init", type=str, default=None, required=False)
parser.add_argument("--no-run", action="store_true")
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

# -------------------------
# TEMPLATE
# -------------------------

@dataclass(frozen=True)
class BenchmarkJob:
    k: int
    s: int
    thread: int
    data: str
    
    @property
    def manifest(self) -> str:
        return os.path.join(self.data, "manifest.txt")
      
    @property
    def name(self) -> str:
        return f"{os.path.basename(self.data)}_{self.k}_{self.s}_{self.thread}"
      
    @property
    def n(self) -> int:
        return int(self.data.rsplit("_", 1)[1])


params: List[BenchmarkJob] = []

with open(args.config, 'r') as f:
    for line in f:
        k, s, thread, data = line.strip().split("\t")
        
        k = int(k)
        s = int(s)
        thread = int(thread)
        
        params.append(BenchmarkJob(k, s, thread, data.rstrip("/")))

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

def estimate_time(k: int, s: int, threads: int, n: int) -> str:
    if n <= 100:
        return "00:30:00"
    elif n <= 500:
        return "02:00:00"
    else:
        return "05:00:00"


def estimate_mem(k: int, s: int, threads: int, n: int):
    if n <= 100:
        return "50G"
    else:
        return "256G"


# -------------------------
# LOAD TEMPLATE
# -------------------------

with open(Path(__file__).parent / "job_template.sh", "r") as f:
    template = f.read()

# -------------------------
# GENERATE JOBS
# -------------------------

job_count = 0
submitted_count = 0

for i, pset in enumerate(params):
    job_id = f"run_{i:04d}"
    job_results = results / job_id

    # Skip completed runs
    launch = not (job_results / "done.flag").exists()
    
    if not launch:
        stat = (job_results / "done.flag").read_text().strip()
        if (stat == "submitted"):
            if args.force:
                launch = True
            else:
                print(f"Skipping {job_id} (already submitted, to override add `--force`)")
                continue
        else:
            assert stat == "complete"
            
            print(f"Skipping {job_id} (already done)")
            continue

    job_results.mkdir(parents=True, exist_ok=True)

    time = estimate_time(pset.k, pset.s, pset.thread, pset.n)
    mem = estimate_mem(pset.k, pset.s, pset.thread, pset.n)

    script_content = template.format(
        job_name=job_id,
        slurm_out=slurm_out,
        threads=pset.thread,
        time=time,
        mem=mem,
        module_load=args.env_init or "",
        outdir=str(job_results),
        manifest=pset.manifest,
        k=pset.k,
        s=pset.s,
        n=pset.n
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
