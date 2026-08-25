from __future__ import annotations

import bz2
import json
import uuid
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence


# ============================================================================
# Clustering Operator Interface
# ============================================================================

class ClusterRefinementOperator(ABC):
    """
    Interface for an object that refines a cluster into one or more
    subclusters.

    Implementations are responsible for determining how the provided
    accessions should be partitioned.
    """
    
    def __init__(self, tag: str) -> None:
        super().__init__()
        
        self.tag = tag

    @abstractmethod
    def refine(self, accessions: Sequence[str]) -> List[List[str]]:
        """
        Split a cluster into subclusters.

        Parameters
        ----------
        accessions
            Members of the cluster being refined.

        Returns
        -------
        List[List[str]]
            A list of subclusters. Each subcluster is a list of accession
            strings.

        Notes
        -----
        - Every accession should appear exactly once in the result.
        - Returning a single cluster identical to the input is treated as
          "no change".
        """
        raise NotImplementedError


class DictClusterRefinement(ClusterRefinementOperator):
    def __init__(self, tag: str, data: dict[str, str]) -> None:
        super().__init__(tag)
        
        self.data = data
    
    def refine(self, accessions: Sequence[str]):
        # Check accessions are known
        missing_accessions = {accn for accn in accessions if accn not in self.data}
        
        if missing_accessions:
            raise ValueError(f"Unknown accessions: {",".join(sorted(missing_accessions))}")
        
        # Group using dictionary
        grouped: dict[str, list[str]] = {}
        
        for accn in accessions:
            grouped.setdefault(self.data[accn], []).append(accn)
        
        return list(grouped.values())


# ============================================================================
# Internal Data Model
# ============================================================================

@dataclass
class ClusterNode:
    cluster_id: str
    tag: str
    accessions: List[str]
    parent_id: Optional[str] = None
    child_ids: List[str] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.child_ids


# ============================================================================
# Cluster Manager
# ============================================================================

class ClusterManager:
    """
    Manages a hierarchical clustering structure.

    The current clustering state is defined by all leaf nodes.
    Non-leaf nodes represent previous clustering levels.
    """

    def __init__(
        self,
        storage_path: str | Path,
        accessions: Optional[Sequence[str]] = None,
    ):
        self.storage_path = Path(storage_path)

        if self.storage_path.exists():
            self._load()
        else:
            if accessions is None:
                raise ValueError(
                    "accessions must be supplied when creating a new cluster file"
                )

            self._create_new("base", accessions)
            self.save()

    # ---------------------------------------------------------------------
    # Persistence
    # ---------------------------------------------------------------------

    def save(self) -> None:
        data = {
            "root_id": self.root_id,
            "nodes": {
                node_id: asdict(node)
                for node_id, node in self.nodes.items()
            },
        }

        with bz2.open(self.storage_path, "wt", encoding="utf-8") as fh:
            json.dump(data, fh, indent=2)

    def _load(self) -> None:
        with bz2.open(self.storage_path, "rt", encoding="utf-8") as fh:
            data = json.load(fh)

        self.root_id = data["root_id"]

        self.nodes: Dict[str, ClusterNode] = {
            node_id: ClusterNode(**node_data)
            for node_id, node_data in data["nodes"].items()
        }

    def _create_new(self, tag: str, accessions: Sequence[str]) -> None:
        root_id = self._new_cluster_id()

        self.nodes: Dict[str, ClusterNode] = {
            root_id: ClusterNode(
                cluster_id=root_id,
                tag=tag,
                accessions=sorted(set(accessions)),
            )
        }

        self.root_id = root_id

    # ---------------------------------------------------------------------
    # Query Methods
    # ---------------------------------------------------------------------

    def iter_clusters(self) -> Iterator[List[str]]:
        """
        Iterate over the current clustering.

        Current clusters are the leaf nodes in the hierarchy.
        """
        for node in self._leaf_nodes():
            yield list(node.accessions)

    def iter_cluster_nodes(self):
        """
        Iterate over current leaf cluster nodes.
        """
        yield from self._leaf_nodes()
        
    def iter_tag(self, tag: str):
        yield from filter(lambda cl: cl.tag == tag, self.nodes.values())

    def get_cluster(self, cluster_id: str) -> List[str]:
        return list(self.nodes[cluster_id].accessions)

    def cluster_count(self) -> int:
        return sum(1 for _ in self._leaf_nodes())

    # ---------------------------------------------------------------------
    # Refinement
    # ---------------------------------------------------------------------

    def refine_clusters(
        self,
        operator: ClusterRefinementOperator,
        persist: bool = True,
    ) -> None:
        """
        Refine every current cluster using the provided operator.

        Each leaf cluster is independently refined. If refinement produces
        multiple subclusters, children are added beneath the cluster node.
        """

        leaves = list(self._leaf_nodes())

        for leaf in leaves:
            self._refine_node(leaf, operator)

        if persist:
            self.save()

    def refine_cluster(
        self,
        cluster_id: str,
        operator: ClusterRefinementOperator,
        persist: bool = True,
    ) -> None:
        """
        Refine a single cluster.
        """

        node = self.nodes[cluster_id]

        if not node.is_leaf:
            raise ValueError(
                f"Cluster '{cluster_id}' has already been refined."
            )

        self._refine_node(node, operator)

        if persist:
            self.save()

    # ---------------------------------------------------------------------
    # Internal
    # ---------------------------------------------------------------------

    def _refine_node(
        self,
        node: ClusterNode,
        operator: ClusterRefinementOperator,
    ) -> None:
        subclusters = operator.refine(node.accessions)

        self._validate_refinement(node.accessions, subclusters)

        # No-op refinement
        if len(subclusters) == 1 and set(subclusters[0]) == set(node.accessions):
            return

        for members in subclusters:
            child_id = self._new_cluster_id()

            child = ClusterNode(
                cluster_id=child_id,
                tag=operator.tag,
                accessions=list(members),
                parent_id=node.cluster_id,
            )

            self.nodes[child_id] = child
            node.child_ids.append(child_id)

    def _leaf_nodes(self) -> Iterator[ClusterNode]:
        for node in self.nodes.values():
            if node.is_leaf:
                yield node

    @staticmethod
    def _validate_refinement(
        original: Sequence[str],
        subclusters: Sequence[Sequence[str]],
    ) -> None:
        original_set = set(original)

        flattened = []
        for cluster in subclusters:
            flattened.extend(cluster)

        refined_set = set(flattened)

        if refined_set != original_set:
            raise ValueError(
                "Refinement must contain exactly the same accessions as "
                "the original cluster."
            )

        if len(flattened) != len(original):
            raise ValueError(
                "Refinement contains duplicates or missing accessions."
            )

    @staticmethod
    def _new_cluster_id() -> str:
        return str(uuid.uuid4())