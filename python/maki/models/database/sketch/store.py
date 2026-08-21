from __future__ import annotations

from dataclasses import asdict
import hashlib
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import json
import os
from pathlib import Path
from typing import Optional

from maki.utils.io import read_fasta, read_fasta_from_gff
from maki.models.database.manifest import DatabasePackage, GenomeData
from sourmash import SourmashSignature

from .core import SketchParameters, load_one_signature, new_minhash, save_one_signature, validate_signature


def _compute_storage_path(accession: str) -> tuple[str, str]:
    """Create a deterministic sharded path with .bz2 extension.
    """
    h = hashlib.sha256(accession.encode()).hexdigest()

    shard1 = h[:2]
    shard2 = h[2:4]

    relpath = f"{shard1}/{shard2}/{accession}.sig.bz2"

    return h, relpath


def _worker_create_signature(
    data: tuple[GenomeData, SketchParameters, str | Path],
) -> tuple[str, str]:
    """
    Worker process.
    """
    item, params, storage_root = data
    
    root = Path(storage_root)

    _, relpath = _compute_storage_path(item.accession)

    outfile = root / "sigs" / relpath
    outfile.parent.mkdir(parents=True, exist_ok=True)
    
    # setch file
    
    mh = new_minhash(params)
    if item.fasta:
        fname = item.fasta
        it = read_fasta(item.fasta)
    elif item.gff:
        fname = item.gff
        it = read_fasta_from_gff(item.gff)
    else:
        raise RuntimeError(f"Genome record {item.accession} did not supply FASTA or GFF data.")
    
    for _, sequence in it:
        # force=True skips k-mers containing ambiguous characters rather than
        # failing the complete genome.
        mh.add_sequence(sequence, force=True)

    if len(mh) == 0:
        raise ValueError(
            f"Genome {item.accession!r} produced an empty sketch"
        )

    signature = SourmashSignature(
        mh,
        name=item.accession,
        filename=str(fname),
    )

    save_one_signature(signature, outfile)
    
    return item.accession, str(relpath)


class SourmashSketchStore:
    """
    Large-scale sourmash signature manager.

    Features
    --------
    * SHA256 sharded filesystem layout
    * gzip-compressed signatures
    * SQLite index
    * multiprocessing batch sketch generation
    * random signature loading
    """

    def __init__(
        self,
        root_dir: str | Path,
        *,
        params: Optional[SketchParameters] = None
    ):
        self.root = Path(root_dir)
        self.root.mkdir(parents=True, exist_ok=True)
        
        self.sigroot = self.root / "sigs"
        self.dbfile = self.root / "accessions.csv"

        self.params = self._read_sketch_parameters(params)
        self.records = self._read_manifest()

        self.params.validate()
        
        self.sigroot.mkdir(parents=True, exist_ok=True)
    
    def _read_sketch_parameters(self, params: Optional[SketchParameters]):
        param_file = self.root / "params.json"
        
        if not param_file.exists():
            if not params:
                raise RuntimeError(f"Missing parameter file at {param_file}")
            else:
                with param_file.open("wt") as f:
                    json.dump(asdict(params), f)
                
                return params
        
        with param_file.open("rt") as f:
            current_params = SketchParameters(**json.load(f))
            
        if params and (current_params != params):
            raise RuntimeError(f"Parameter mismatch! Existing params: {asdict(current_params)}")
            
        return current_params
        
    def _read_manifest(self):
        recs: dict[str, str] = {}
        
        if self.dbfile.exists():
            with self.dbfile.open("rt") as f:
                for line in f:
                    data = line.strip().split(",")
                    if data:
                        recs[data[0]] = data[1]
                    
        return recs
    
    def _write_manifest(self):
        with self.dbfile.open("wt") as f:
            for accn, sig_path in self.records.items():
                f.write(f"{accn},{sig_path}\n")

    def accessions(self):
        yield from self.records

    # --------------------------------------------------
    # batch sketching
    # --------------------------------------------------

    def sketch_many(
        self,
        package: DatabasePackage,
        *,
        workers: int | None = None,
        skip_existing: bool = False
    ):
        """Sketch many FASTA files in parallel.
        """
        # Check existing records
        
        new_accessions: set[str] = set()
        existing_accessions: set[str] = set()
        
        for rec in package.records():
            if rec.accession in self.records:
                existing_accessions.add(rec.accession)
            else:
                new_accessions.add(rec.accession)
        
        if existing_accessions and not skip_existing:
            raise RuntimeError(f"{len(existing_accessions)} accessions already in sketch store, use `--skip-existing` if you want to skip these.")
                
        # Import sketches
        
        with package.retrieve_data(new_accessions) as data:
            sketch_jobs = zip(data.genomes(), repeat(self.params), repeat(self.root))
            
            with ProcessPoolExecutor(max_workers=workers) as executor:
                sketched = list(executor.map(_worker_create_signature, sketch_jobs))
                
        # Update manifest
        
        for accn, sig_file in sketched:
            self.records[accn] = sig_file
            
        self._write_manifest()

    # --------------------------------------------------
    # loading
    # --------------------------------------------------

    def load_signature(
        self,
        accession: str,
    ) -> SourmashSignature:
        """
        Load a SourmashSignature object.
        """
        sketch_path = self.sigroot / self.records[accession]
        
        sig = load_one_signature(sketch_path)
        validate_signature(sig, self.params, allow_finer_scaled=False)
        
        return sig

    def load_many(
        self,
        accessions: list[str],
    ) -> dict[str, SourmashSignature]:
        return {
            acc: self.load_signature(acc)
            for acc in accessions
        }
