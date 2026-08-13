from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, List, Sequence, Tuple


class BenchmarkPhase(ABC):
    """
    A single measurable stage within a benchmark workflow.

    The framework measures execution time and memory.
    Implementations only perform computation.
    """

    name: str

    @abstractmethod
    def transform(
        self,
        data: Any,
    ) -> Tuple[List[str], Any]:
        """
        Execute the phase.

        Parameters
        ----------
        data
            Input object produced by the previous phase.

        Returns
        -------
        Any
            Output passed to the next phase.
        """
        

class BenchmarkWorkflow(ABC):
    """
    Ordered sequence of benchmark phases.

    The framework executes phases sequentially
    and automatically propagates outputs.
    """
    
    registry = {}
    
    name: str
    phase_names: list[str]
    compatible_data_types: list[str]
    
    def __init_subclass__(cls, **kwargs) -> None:
        super().__init_subclass__(**kwargs)

        name = getattr(cls, "name", None)
        if isinstance(name, str):
            BenchmarkWorkflow.registry[name] = cls
    
    @classmethod
    def from_name(cls, name: str, workflow_args: dict) -> "BenchmarkWorkflow":
        try:
            workflow_type = BenchmarkWorkflow.registry[name]
        except KeyError:
            raise TypeError(f"Unrecognosed workflow name: {name}")
        
        workflow_args, missing_args, extra_args = workflow_type.validate_args(workflow_args)
        if missing_args:
            raise RuntimeError(f"Workflow {name} requires the following missing arguments: {", ".join(f"{arg}<{arg_type}>" for arg, arg_type in missing_args)}")
        elif extra_args:
            raise RuntimeError(f"Workflow {name} received unrecognised arguments: {", ".join(extra_args)}")
        else:
            return workflow_type(**workflow_args)
    
    @staticmethod
    @abstractmethod
    def validate_args(args: dict) -> Tuple[dict, List[Tuple[str, str]], List[str]]:
        """
        Identify missing CLI arguments.
        """

    @abstractmethod
    def phases(self) -> Sequence[BenchmarkPhase]:
        """
        Return ordered workflow phases.
        """
        