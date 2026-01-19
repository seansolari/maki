
export NEXTFLOW_PATH=nextflow
export MAKI_PATH=/mnt/BSTORE/ssol0002/SOFTWARE/pam2s/build/bin/maki
export MAKI_PIPELINE=/mnt/BSTORE/ssol0002/SOFTWARE/pam2s/maki_nf/main.nf
export MAKI_TEST_PROFILE=conda      # or docker/standard
export MAKI_TEST_TIMEOUT=30         # override if needed
export MAKI_TEST_SEED=42            # to make results reproducible

pytest -q
