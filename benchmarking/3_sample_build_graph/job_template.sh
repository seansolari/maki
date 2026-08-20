#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={slurm_out}/{job_name}.out
#SBATCH --error={slurm_out}/{job_name}.err
{email_region}
#SBATCH --tasks=1
#SBATCH --cpus-per-task={threads}
#SBATCH --time=04:00:00
#SBATCH --mem=64G

set -euo pipefail

echo "Running job: {job_name}"
echo "Threads: {threads}"

mkdir -p {outdir}

{module_load}

# Run benchmark

maki benchmark run-pe \
    --forward {forward} --reverse {reverse} {subsample_region} \
    --stats-out {outdir}/stats.json \
    --workflow graph_build \
    -k {k} -s {s} --threads {threads} --graph-out {outdir}/index \
    > {outdir}/stdout.txt
    2> {outdir}/stderr.txt

# Size of database on disk
du -bsh {outdir}/index > {outdir}/index-space.txt

# Remove index
rm -rf {outdir}/index

# Diagnostics
echo {forward} {reverse} {subsample} {k} {s} {threads} > {outdir}/params.txt
hostname > {outdir}/host.txt
date > {outdir}/date.txt

# Mark completion
echo 'complete' > {outdir}/done.flag
