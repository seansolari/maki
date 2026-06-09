
import gzip
import json
import logging
import os
from pathlib import Path
import subprocess
from typing import List, Sequence, Tuple, Union

import bakta.io.gff as gff
import bakta.config as cfg
from xopen import xopen
from tqdm.contrib.concurrent import process_map


logger = logging.getLogger(__name__)


def run_and_gzip_output(
    cmd: Union[str, Sequence[str]],
    stdout_gz_path: Path,
    stderr_gz_path: Path
) -> None:
    """
    Run a subprocess command and write stdout and stderr to gzip files.

    Args:
        cmd: Command to execute (list of args or string if shell=True)
        stdout_gz_path: Path to gzip file for stdout
        stderr_gz_path: Path to gzip file for stderr
    Raises:
        subprocess.CalledProcessError: if the command exits non-zero
    """

    # Open gzip files in binary mode
    with gzip.open(stdout_gz_path, "wb") as gz_out, \
         gzip.open(stderr_gz_path, "wb") as gz_err:

        # Start the process
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        assert process.stdout
        assert process.stderr

        # Stream output incrementally (avoids buffering everything in memory)
        try:
            while True:
                out = process.stdout.readline()
                err = process.stderr.readline()

                if out:
                    gz_out.write(out)
                if err:
                    gz_err.write(err)

                # Exit when process is done and no more output
                if not out and not err and process.poll() is not None:
                    break

            returncode = process.wait()

        finally:
            if process.stdout:
                process.stdout.close()
            if process.stderr:
                process.stderr.close()

        # Mimic subprocess.check_call behavior
        if returncode != 0:
            raise subprocess.CalledProcessError(returncode, cmd)
          

def bakta_to_json(json_path: Path, gff_path: Path):
    with xopen(json_path, threads=0) as fh:
        data = json.load(fh)
    
    cfg.db_info = {
        'major': data['version']['db']['version'].split('.')[0],
        'minor': data['version']['db']['version'].split('.')[1],
        'type': data['version']['db']['type']
    }
    
    features_by_sequence = {seq['id']: [] for seq in data['sequences']}
    for feature in data['features']:
        seq_id = feature['sequence'] if 'sequence' in feature else feature['contig']  # <1.10.0 compatibility
        features_by_sequence[seq_id].append(feature)

    gff.write_features(data, features_by_sequence, gff_path) # type: ignore


def json_to_gff(json: Path, gff: Path) -> None:
    if not gff.exists():
        gff.parent.mkdir(parents=True, exist_ok=True)
        
        gff_temp = gff.with_suffix(".gff.temp")
        
        try:
            bakta_to_json(json, gff_temp)
            os.replace(gff_temp, gff)
        finally:
            gff_temp.unlink(missing_ok=True)


def pjson_to_gff(jobs: List[Tuple[Path, Path]], concurrency: int):    
    process_map(json_to_gff, *zip(*jobs), max_workers=concurrency, chunksize=1, desc="Writing GFF")
