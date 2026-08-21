"""
three_pass_sourmash_cluster.py

Online representative-based clustering for scaled sourmash MinHash sketches.

Strategy
--------
Pass 1:
    Compare the query against the union of all representative index sketches.
    If query-in-union containment is below screen_containment_threshold,
    create a new cluster immediately.

Pass 2:
    Look up each coarse/downsampled query hash in an inverted index.
    Score all resulting representative candidates using verification sketches.

Pass 3:
    If no indexed candidate is found, compare exhaustively against all
    representatives. By default, exhaustive comparison also occurs when
    candidates exist but none pass the final threshold, because that is safer.

Requirements
------------
    sourmash >= 4
    Scaled/FracMinHash sketches only; fixed-size ``num`` sketches are rejected.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import (
    Callable,
    DefaultDict,
    Dict,
    Hashable,
    Iterable,
    List,
    Mapping,
    Optional,
    Set,
    Tuple,
    Union,
)

from sourmash import MinHash


MinHashLike = Union[MinHash, object]
RepresentativeID = Hashable
ScoreFunction = Callable[[MinHash, MinHash], float]


@dataclass(frozen=True)
class Assignment:
    """
    Result returned by ThreePassClusterer.assign().

    Attributes
    ----------
    representative_id:
        Identifier of the selected or newly created representative.
    is_new_cluster:
        True when the query became a new representative.
    score:
        Final query/representative similarity. It is 1.0 for a newly
        created representative.
    union_containment:
        Query-in-union containment observed during pass 1.
    candidates_examined:
        Number of representatives compared using full verification sketches.
    index_candidates:
        Number of distinct representatives returned by the inverted index.
    pass_used:
        One of:
            "initial"
            "union-negative"
            "inverted-index"
            "exhaustive-no-index-match"
            "exhaustive-candidate-failure"
            "new-after-index"
            "new-after-exhaustive"
    """

    representative_id: RepresentativeID
    is_new_cluster: bool
    score: float
    union_containment: float
    candidates_examined: int
    index_candidates: int
    pass_used: str


def get_minhash(value: MinHashLike) -> MinHash:
    """
    Return a MinHash from either a MinHash or SourmashSignature-like object.

    A SourmashSignature exposes its sketch through ``.minhash``.
    """
    if isinstance(value, MinHash):
        return value

    mh = getattr(value, "minhash", None)
    if mh is None:
        raise TypeError(
            "Expected a sourmash.MinHash or an object exposing '.minhash'."
        )
    return mh


def flatten_minhash(mh: MinHash) -> MinHash:
    """
    Return an abundance-free sketch suitable for set-based comparisons.
    """
    if getattr(mh, "track_abundance", False):
        return mh.flatten()
    return mh


def downsample_scaled(mh: MinHash, scaled: int) -> MinHash:
    """
    Flatten and downsample a scaled MinHash to ``scaled``.

    A larger scaled value is a coarser sketch. Upsampling is impossible,
    so ``scaled`` must be at least the sketch's current scaled value.
    """
    mh = flatten_minhash(mh)

    if not mh.scaled:
        raise ValueError(
            "This implementation requires scaled/FracMinHash sketches; "
            "fixed-size num sketches are not supported."
        )

    if scaled < mh.scaled:
        raise ValueError(
            f"Cannot upsample from scaled={mh.scaled} to scaled={scaled}. "
            "Choose a target scaled value greater than or equal to every "
            "input sketch's scaled value."
        )

    if scaled == mh.scaled:
        return mh

    return mh.downsample(scaled=scaled)


def hash_values(mh: MinHash) -> Set:
    """
    Return the retained hash values from a MinHash as a Python set.
    """
    return set(mh.hashes)


def jaccard_score(query: MinHash, representative: MinHash) -> float:
    """
    Compute set Jaccard between already-compatible sketches.
    """
    q = hash_values(query)
    r = hash_values(representative)

    union_size = len(q) + len(r) - len(q & r)
    if union_size == 0:
        return 1.0

    return len(q & r) / union_size


def query_containment_score(
    query: MinHash,
    representative: MinHash,
) -> float:
    """
    Compute query-in-representative containment:

        |Q intersect R| / |Q|
    """
    q = hash_values(query)
    if not q:
        return 0.0

    r = hash_values(representative)
    return len(q & r) / len(q)


def max_containment_score(
    query: MinHash,
    representative: MinHash,
) -> float:
    """
    Compute maximum bidirectional containment.
    """
    q = hash_values(query)
    r = hash_values(representative)
    intersection = len(q & r)

    q_containment = intersection / len(q) if q else 0.0
    r_containment = intersection / len(r) if r else 0.0

    return max(q_containment, r_containment)


def average_containment_score(
    query: MinHash,
    representative: MinHash,
) -> float:
    """
    Compute the mean of the two directional containments.
    """
    q = hash_values(query)
    r = hash_values(representative)
    intersection = len(q & r)

    q_containment = intersection / len(q) if q else 0.0
    r_containment = intersection / len(r) if r else 0.0

    return 0.5 * (q_containment + r_containment)


def make_sourmash_ani_scorer(
    *,
    kind: str = "average",
    prob_threshold: float = 0.001,
) -> ScoreFunction:
    """
    Construct a scorer using sourmash ANI estimation methods.

    Parameters
    ----------
    kind:
        "average":
            Use average-containment ANI.
        "max":
            Use maximum-containment ANI.

    prob_threshold:
        Estimation confidence/probability cutoff passed to sourmash.

    Notes
    -----
    Sourmash versions may return either a float-like value or a result
    object carrying an ``ani`` attribute. This wrapper handles both.
    """

    def extract_ani(result: object) -> float:
        if hasattr(result, "ani"):
            value = getattr(result, "ani")
        else:
            value = result

        if value is None:
            return 0.0

        try:
            value = float(value)
        except (TypeError, ValueError):
            return 0.0

        # Treat NaN as a failed estimate.
        if value != value:
            return 0.0

        return value

    if kind == "average":

        def scorer(query: MinHash, representative: MinHash) -> float:
            result = query.avg_containment_ani(
                representative,
                downsample=False,
                prob_threshold=prob_threshold,
            )
            return extract_ani(result)

        return scorer

    if kind == "max":

        def scorer(query: MinHash, representative: MinHash) -> float:
            result = query.max_containment_ani(
                representative,
                downsample=False,
                prob_threshold=prob_threshold,
            )
            return extract_ani(result)

        return scorer

    raise ValueError("kind must be either 'average' or 'max'.")


class ThreePassClusterer:
    """
    Incremental three-pass clusterer for scaled sourmash sketches.

    Parameters
    ----------
    verification_scaled:
        Scaled value used for final representative comparisons.

        This must be greater than or equal to the scaled value of every
        incoming signature. A smaller value would require unavailable hashes.

    index_scaled:
        Coarser scaled value used by the inverted index and union filter.
        It must be >= verification_scaled.

        Larger values reduce memory and posting-list work, but increase the
        probability that a related pair shares no indexed hash.

    screen_containment_threshold:
        Minimum query-in-union containment required to proceed beyond pass 1.

        For a final query-containment threshold, using the same value is safe.
        For Jaccard threshold t, using t is also a safe negative filter because

            query containment >= Jaccard.

        For ANI clustering, do not set this directly to the ANI threshold.
        ANI and containment are on different scales. Use a containment cutoff
        corresponding to the ANI/k-mer/sketch configuration.

    match_threshold:
        Minimum score from ``score_fn`` required for cluster assignment.

    score_fn:
        Function accepting two compatible verification MinHashes and returning
        a score where larger is better. Defaults to set Jaccard.

    exhaustive_on_candidate_failure:
        If True, run exhaustive comparison when indexed candidates exist but
        none passes match_threshold. This protects against false negatives
        caused by the coarse inverted index.

        If False, exhaustive fallback happens only when the inverted index
        returns zero representatives, exactly matching the narrower strategy
        described in the question.

    require_index_overlap:
        Minimum number of coarse indexed hashes shared by a representative
        before it becomes a pass-2 candidate. The default of 1 maximizes
        sensitivity.
    """

    def __init__(
        self,
        *,
        verification_scaled: int,
        index_scaled: int,
        screen_containment_threshold: float,
        match_threshold: float,
        score_fn: ScoreFunction = jaccard_score,
        exhaustive_on_candidate_failure: bool = True,
        require_index_overlap: int = 1,
    ) -> None:
        if verification_scaled <= 0:
            raise ValueError("verification_scaled must be positive.")

        if index_scaled < verification_scaled:
            raise ValueError(
                "index_scaled must be greater than or equal to "
                "verification_scaled."
            )

        if not 0.0 <= screen_containment_threshold <= 1.0:
            raise ValueError(
                "screen_containment_threshold must be between 0 and 1."
            )

        if not 0.0 <= match_threshold <= 1.0:
            raise ValueError("match_threshold must be between 0 and 1.")

        if require_index_overlap < 1:
            raise ValueError("require_index_overlap must be at least 1.")

        self.verification_scaled = verification_scaled
        self.index_scaled = index_scaled
        self.screen_containment_threshold = (
            screen_containment_threshold
        )
        self.match_threshold = match_threshold
        self.score_fn = score_fn
        self.exhaustive_on_candidate_failure = (
            exhaustive_on_candidate_failure
        )
        self.require_index_overlap = require_index_overlap

        # Full verification sketches, normalized to verification_scaled.
        self.representatives: Dict[RepresentativeID, MinHash] = {}

        # Coarse index sketches. Keeping these makes removal/checkpointing easy.
        self.index_sketches: Dict[RepresentativeID, Set[int]] = {}

        # Coarse hash -> IDs of representatives containing that hash.
        self.inverted_index: DefaultDict[
            int, Set[RepresentativeID]
        ] = defaultdict(set)

        # Union of all coarse representative hashes.
        self.representative_union: Set[int] = set()

        # Reference sketch characteristics established by the first input.
        self._compatibility_key: Optional[Tuple[object, ...]] = None

        # Used only when the caller does not provide representative IDs.
        self._next_integer_id = 0

    def __len__(self) -> int:
        return len(self.representatives)

    @staticmethod
    def _sketch_compatibility_key(mh: MinHash) -> Tuple[object, ...]:
        """
        Identify sketch characteristics that must agree across comparisons.
        """
        return (
            mh.ksize,
            getattr(mh, "seed", None),
            getattr(mh, "moltype", None),
            getattr(mh, "is_protein", None),
            getattr(mh, "dayhoff", None),
            getattr(mh, "hp", None),
            getattr(mh, "skipm1n3", None),
            getattr(mh, "skipm2n3", None),
        )

    def _normalize_query(
        self,
        value: MinHashLike,
    ) -> Tuple[MinHash, Set[int]]:
        """
        Return verification sketch and coarse-index hashes for an input.
        """
        original = get_minhash(value)
        verification = downsample_scaled(
            original,
            self.verification_scaled,
        )

        compatibility_key = self._sketch_compatibility_key(verification)
        if self._compatibility_key is None:
            self._compatibility_key = compatibility_key
        elif compatibility_key != self._compatibility_key:
            raise ValueError(
                "Input sketch is incompatible with existing representatives. "
                "Check k-mer size, molecule type, hash seed, and skip-mer/"
                "protein encoding."
            )

        indexed = downsample_scaled(verification, self.index_scaled)
        return verification, hash_values(indexed)

    def _allocate_id(self) -> int:
        while self._next_integer_id in self.representatives:
            self._next_integer_id += 1

        result = self._next_integer_id
        self._next_integer_id += 1
        return result

    def add_representative(
        self,
        signature: MinHashLike,
        representative_id: Optional[RepresentativeID] = None,
    ) -> RepresentativeID:
        """
        Add a representative directly to the clusterer.

        This updates:
            * the verification-sketch mapping;
            * the coarse inverted index;
            * the accumulated representative-union sketch.
        """
        verification, index_hashes = self._normalize_query(signature)

        if representative_id is None:
            representative_id = self._allocate_id()

        if representative_id in self.representatives:
            raise KeyError(
                f"Representative ID already exists: {representative_id!r}"
            )

        self.representatives[representative_id] = verification
        self.index_sketches[representative_id] = index_hashes

        for hash_value in index_hashes:
            self.inverted_index[hash_value].add(representative_id)

        self.representative_union.update(index_hashes)
        return representative_id

    def union_containment_from_hashes(
        self,
        query_index_hashes: Set[int],
    ) -> float:
        """
        Compute query-in-union containment at index_scaled.
        """
        if not query_index_hashes:
            return 0.0

        common = len(
            query_index_hashes.intersection(self.representative_union)
        )
        return common / len(query_index_hashes)

    def find_index_candidates(
        self,
        query_index_hashes: Set[int],
    ) -> Dict[RepresentativeID, int]:
        """
        Return representatives sharing indexed hashes with the query.

        The returned mapping is:

            representative_id -> number of shared index hashes
        """
        overlap_counts: DefaultDict[RepresentativeID, int]
        overlap_counts = defaultdict(int)

        for hash_value in query_index_hashes:
            for representative_id in self.inverted_index.get(
                hash_value, ()
            ):
                overlap_counts[representative_id] += 1

        return {
            representative_id: count
            for representative_id, count in overlap_counts.items()
            if count >= self.require_index_overlap
        }

    def best_match(
        self,
        query: MinHash,
        candidate_ids: Iterable[RepresentativeID],
    ) -> Tuple[Optional[RepresentativeID], float, int]:
        """
        Return the highest-scoring candidate.

        Ties are resolved deterministically by retaining the first candidate
        encountered. Candidate IDs are sorted by repr() to avoid relying on
        cross-type comparison between arbitrary hashable identifiers.
        """
        unique_ids = set(candidate_ids)
        ordered_ids = sorted(unique_ids, key=repr)

        best_id: Optional[RepresentativeID] = None
        best_score = float("-inf")
        examined = 0

        for representative_id in ordered_ids:
            representative = self.representatives[representative_id]
            score = float(self.score_fn(query, representative))
            examined += 1

            if score > best_score:
                best_id = representative_id
                best_score = score

        if best_id is None:
            return None, 0.0, examined

        return best_id, best_score, examined

    def _create_cluster(
        self,
        verification_query: MinHash,
        *,
        representative_id: Optional[RepresentativeID],
        union_containment: float,
        candidates_examined: int,
        index_candidates: int,
        pass_used: str,
    ) -> Assignment:
        new_id = self.add_representative(
            verification_query,
            representative_id=representative_id,
        )

        return Assignment(
            representative_id=new_id,
            is_new_cluster=True,
            score=1.0,
            union_containment=union_containment,
            candidates_examined=candidates_examined,
            index_candidates=index_candidates,
            pass_used=pass_used,
        )

    def assign(
        self,
        signature: MinHashLike,
        *,
        new_representative_id: Optional[RepresentativeID] = None,
    ) -> Assignment:
        """
        Assign one signature to an existing representative or create a cluster.

        Existing representatives are not modified when a query joins a cluster.
        This is therefore leader clustering, not centroid updating.
        """
        verification_query, query_index_hashes = self._normalize_query(
            signature
        )

        # First object always creates the initial cluster.
        if not self.representatives:
            return self._create_cluster(
                verification_query,
                representative_id=new_representative_id,
                union_containment=0.0,
                candidates_examined=0,
                index_candidates=0,
                pass_used="initial",
            )

        # Pass 1: query containment in the union of all representatives.
        union_containment = self.union_containment_from_hashes(
            query_index_hashes
        )

        if union_containment < self.screen_containment_threshold:
            return self._create_cluster(
                verification_query,
                representative_id=new_representative_id,
                union_containment=union_containment,
                candidates_examined=0,
                index_candidates=0,
                pass_used="union-negative",
            )

        # Pass 2: retrieve representatives through the inverted index.
        overlap_counts = self.find_index_candidates(query_index_hashes)
        candidate_ids = set(overlap_counts)
        n_index_candidates = len(candidate_ids)

        if candidate_ids:
            best_id, best_score, examined = self.best_match(
                verification_query,
                candidate_ids,
            )

            if best_id is not None and best_score >= self.match_threshold:
                return Assignment(
                    representative_id=best_id,
                    is_new_cluster=False,
                    score=best_score,
                    union_containment=union_containment,
                    candidates_examined=examined,
                    index_candidates=n_index_candidates,
                    pass_used="inverted-index",
                )

            if not self.exhaustive_on_candidate_failure:
                return self._create_cluster(
                    verification_query,
                    representative_id=new_representative_id,
                    union_containment=union_containment,
                    candidates_examined=examined,
                    index_candidates=n_index_candidates,
                    pass_used="new-after-index",
                )

            # Safer fallback: compare only representatives not already checked.
            remaining_ids = (
                set(self.representatives).difference(candidate_ids)
            )
            fallback_id, fallback_score, fallback_examined = self.best_match(
                verification_query,
                remaining_ids,
            )
            total_examined = examined + fallback_examined

            # The best indexed score may still exceed every fallback score,
            # even though it did not pass the assignment threshold.
            if fallback_id is not None and fallback_score > best_score:
                best_id = fallback_id
                best_score = fallback_score

            if best_id is not None and best_score >= self.match_threshold:
                return Assignment(
                    representative_id=best_id,
                    is_new_cluster=False,
                    score=best_score,
                    union_containment=union_containment,
                    candidates_examined=total_examined,
                    index_candidates=n_index_candidates,
                    pass_used="exhaustive-candidate-failure",
                )

            return self._create_cluster(
                verification_query,
                representative_id=new_representative_id,
                union_containment=union_containment,
                candidates_examined=total_examined,
                index_candidates=n_index_candidates,
                pass_used="new-after-exhaustive",
            )

        # Pass 3: no representative shared an indexed hash.
        best_id, best_score, examined = self.best_match(
            verification_query,
            self.representatives,
        )

        if best_id is not None and best_score >= self.match_threshold:
            return Assignment(
                representative_id=best_id,
                is_new_cluster=False,
                score=best_score,
                union_containment=union_containment,
                candidates_examined=examined,
                index_candidates=0,
                pass_used="exhaustive-no-index-match",
            )

        return self._create_cluster(
            verification_query,
            representative_id=new_representative_id,
            union_containment=union_containment,
            candidates_examined=examined,
            index_candidates=0,
            pass_used="new-after-exhaustive",
        )
