
from dataclasses import asdict, dataclass
import json
import logging
import os
from pathlib import Path
import uuid

from maki.sketch.core import SketchParameters

from .manysearch import ManySearchStore
from .sketch import MetagenomeSketch, MetagenomeSketchStore


logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class ClassifyParameters:
    ...


class OutputHook:
    """Lightweight interface to classification results.
    """
    
    PARAMETERS_FILE = "parameters.json"
    
    # Output preparation and hook initialisation
    # ------------------------------------------

    def __init__(
        self,
        output_root: str | Path,
    ) -> None:
        self.root = Path(output_root)

        self.parameters_file = self.root / self.PARAMETERS_FILE
        
        if not self.parameters_file.is_file():
            raise FileNotFoundError(
                f"Classification parameter file does not exist: {self.parameters_file}"
            )
        
        # Read parameters
        
        with self.parameters_file.open("r", encoding="utf-8") as handle:
            all_params = json.load(handle)
        
        self.params = ClassifyParameters(**all_params["classify"])
        
        # Handle to current sketches
        
        sketch_params = SketchParameters(**all_params["sketch"])
        self.sketches = MetagenomeSketchStore(self.root, params=sketch_params)
        
    @classmethod
    def init(
        cls,
        output_root: str | Path,
        params: ClassifyParameters,
        sketch_params: SketchParameters,
        *,
        exist_ok: bool = False
    ):
        root = Path(output_root)
        
        if root.exists():
            if not root.is_dir():
                raise NotADirectoryError(root)

            if not exist_ok:
                raise FileExistsError(root)

            hook = cls(root)
            
            if hook.params != params:
                raise ValueError(
                    "Existing classification parameters do not match "
                    "the supplied parameters"
                )
            
            if hook.sketches.params != sketch_params:
                raise ValueError(
                    "Existing sketch parameters do not match the "
                    "supplied parameters"
                )
                
            return hook
        
        root.parent.mkdir(parents=True, exist_ok=True)
        
        staging_root = root.with_name(
            f".{root.name}.initializing-{uuid.uuid4().hex}"
        )
        
        try:
            staging_root.mkdir()
            
            # create metadata file
            
            all_params = {
                "classify": asdict(params),
                "sketch": asdict(sketch_params)
            }
            
            param_file = staging_root / cls.PARAMETERS_FILE
            
            with param_file.open("w", encoding="utf-8") as fh:
                json.dump(all_params, fh, indent=2, sort_keys=True)
                fh.write("\n")
                fh.flush()
                os.fsync(fh.fileno())
            
            staging_root.rename(root)
        
        except Exception:
            if staging_root.exists():
                import shutil
                
                shutil.rmtree(staging_root, ignore_errors=True)
            
            raise
        
        return cls(root)
    
    def screen(self, rocksdb: Path, parallel: int):
        mgr = ManySearchStore(self.root, rocksdb, self.sketches.params)
        return mgr.update(self.sketches.list_sketches(), threads=parallel)
