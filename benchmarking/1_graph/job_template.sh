#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output=logs/{job_name}.out
#SBATCH --error=logs/{job_name}.err
#SBATCH --tasks=1
#SBATCH --cpus-per-task={threads}
#SBATCH --time={time}
#SBATCH --mem={mem}

set -euo pipefail

echo "Running job: {job_name}"
echo "Threads: {threads}"

mkdir -p {outdir}

{module_load}

# Run benchmark
/usr/bin/time -v -o {outdir}/time.txt \
    bash -c "
        maki to-graph {manifest} {k} {s} {threads} {outdir}/index \
        > {outdir}/stdout.txt \
        2> {outdir}/stderr.txt
    "

# Size of database on disk
du -bsh {outdir}/index > {outdir}/index-space.txt

# Remove index
rm -rf {outdir}/index

# Diagnostics
echo {n} {k} {s} {threads} > {outdir}/params.txt
hostname > {outdir}/host.txt
date > {outdir}/date.txt

# Mark completion
echo 'complete' > {outdir}/done.flag
