"""Immunarch-inspired repertoire summaries with Python-native outputs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from ._optional import require_dependency
from .bulk import abundance_matrix

DEFAULT_PROPORTION_BINS = {
    "Hyperexpanded": 1e-2,
    "Large": 1e-3,
    "Medium": 1e-4,
    "Small": 1e-5,
    "Rare": 1e-6,
}
DEFAULT_RANK_BINS = (10, 30, 100, 300, 1_000, 10_000, 100_000)


def _validate_abundance(abundance: pd.DataFrame) -> None:
    required = {"sample_id", "clone_id", "count"}
    missing = required.difference(abundance.columns)
    if missing:
        raise KeyError(f"Abundance table is missing columns: {sorted(missing)}")
    counts = pd.to_numeric(abundance["count"], errors="coerce")
    if counts.isna().any() or (counts < 0).any():
        raise ValueError("Counts must be finite and non-negative")


def _with_frequency(abundance: pd.DataFrame) -> pd.DataFrame:
    _validate_abundance(abundance)
    result = abundance.copy()
    result["count"] = pd.to_numeric(result["count"], errors="raise").astype(float)
    totals = result.groupby("sample_id", observed=True)["count"].transform("sum")
    result["frequency"] = np.divide(
        result["count"],
        totals,
        out=np.zeros(len(result), dtype=float),
        where=totals.to_numpy() > 0,
    )
    return result


def coverage_diversity(
    abundance: pd.DataFrame,
    *,
    percentages: Sequence[float] = (50,),
) -> pd.DataFrame:
    """Calculate DXX: clones required to occupy each repertoire percentage."""

    requested = [float(value) for value in percentages]
    if not requested or any(value <= 0 or value > 100 for value in requested):
        raise ValueError("percentages must contain values in (0, 100]")
    records: list[dict[str, object]] = []
    for sample_id, sample in _with_frequency(abundance).groupby(
        "sample_id", observed=True, sort=False
    ):
        frequencies = np.sort(sample["frequency"].to_numpy(dtype=float))[::-1]
        cumulative = np.cumsum(frequencies)
        richness = int((frequencies > 0).sum())
        for percentage in requested:
            target = percentage / 100
            dxx = (
                int(np.searchsorted(cumulative + 1e-12, target, side="left") + 1)
                if richness
                else 0
            )
            dxx = min(dxx, richness)
            records.append(
                {
                    "sample_id": str(sample_id),
                    "percentage": percentage,
                    "dxx": dxx,
                    "richness": richness,
                    "dxx_fraction": dxx / richness if richness else np.nan,
                }
            )
    return pd.DataFrame.from_records(records)


def hill_diversity(
    abundance: pd.DataFrame,
    *,
    orders: Sequence[float] = (0, 1, 2, 3, 4, 5),
    round_counts: bool = True,
) -> pd.DataFrame:
    """Calculate a Hill-number diversity profile through scikit-bio."""

    q_values = [float(order) for order in orders]
    if not q_values or any(order < 0 for order in q_values):
        raise ValueError("orders must contain non-negative values")
    diversity = require_dependency(
        "skbio.diversity", extra="bulk", feature="Hill diversity"
    )
    matrix = abundance_matrix(abundance)
    if round_counts:
        matrix = matrix.round()
    records: list[dict[str, object]] = []
    for order in q_values:
        values = diversity.alpha_diversity(
            "hill",
            matrix.to_numpy(),
            ids=matrix.index.astype(str).tolist(),
            order=order,
        )
        records.extend(
            {
                "sample_id": str(sample_id),
                "q": order,
                "hill_number": float(value),
            }
            for sample_id, value in values.items()
        )
    return pd.DataFrame.from_records(records)


def rank_abundance(abundance: pd.DataFrame, *, limit: int | None = None) -> pd.DataFrame:
    """Return within-sample clone ranks and normalized abundance."""

    result = _with_frequency(abundance).sort_values(
        ["sample_id", "frequency", "clone_id"],
        ascending=[True, False, True],
    )
    result["rank"] = result.groupby("sample_id", observed=True).cumcount() + 1
    if limit is not None:
        if limit < 1:
            raise ValueError("limit must be positive")
        result = result[result["rank"] <= limit]
    return result.reset_index(drop=True)


def annotate_clonality_proportion(
    abundance: pd.DataFrame,
    *,
    bins: Mapping[str, float] = DEFAULT_PROPORTION_BINS,
    key_added: str = "clonal_prop_bin",
) -> pd.DataFrame:
    """Annotate clones using immunarch-style abundance classes."""

    thresholds = sorted(
        ((str(label), float(value)) for label, value in bins.items()),
        key=lambda item: item[1],
        reverse=True,
    )
    if not thresholds or any(value <= 0 or value > 1 for _, value in thresholds):
        raise ValueError("bin thresholds must be in (0, 1]")
    result = _with_frequency(abundance)

    def classify(frequency: float) -> str:
        for label, threshold in thresholds:
            if frequency >= threshold:
                return label
        return "Ultra-rare"

    categories = [label for label, _ in thresholds] + ["Ultra-rare"]
    result[key_added] = pd.Categorical(
        result["frequency"].map(classify), categories=categories, ordered=True
    )
    return result


def clonality_proportion(
    abundance: pd.DataFrame,
    *,
    bins: Mapping[str, float] = DEFAULT_PROPORTION_BINS,
) -> pd.DataFrame:
    """Aggregate occupied repertoire space by abundance classes."""

    annotated = annotate_clonality_proportion(abundance, bins=bins)
    return (
        annotated.groupby(
            ["sample_id", "clonal_prop_bin"], observed=False, sort=False
        )
        .agg(
            n_clonotypes=("clone_id", "nunique"),
            occupied_frequency=("frequency", "sum"),
            count=("count", "sum"),
        )
        .reset_index()
    )


def annotate_clonality_rank(
    abundance: pd.DataFrame,
    *,
    bins: Sequence[int] = DEFAULT_RANK_BINS,
    key_added: str = "clonal_rank_bin",
) -> pd.DataFrame:
    """Annotate clones by within-repertoire rank bins."""

    thresholds = sorted({int(value) for value in bins})
    if not thresholds or thresholds[0] < 1:
        raise ValueError("rank bins must be positive integers")
    result = rank_abundance(abundance)
    labels: list[str] = []
    lower = 1
    for threshold in thresholds:
        labels.append(f"{lower}-{threshold}")
        lower = threshold + 1
    labels.append(f">{thresholds[-1]}")
    indices = np.searchsorted(thresholds, result["rank"].to_numpy(), side="left")
    result[key_added] = pd.Categorical(
        [labels[min(index, len(labels) - 1)] for index in indices],
        categories=labels,
        ordered=True,
    )
    return result


def clonality_rank(
    abundance: pd.DataFrame,
    *,
    bins: Sequence[int] = DEFAULT_RANK_BINS,
) -> pd.DataFrame:
    """Aggregate occupied repertoire space by clone-rank bins."""

    annotated = annotate_clonality_rank(abundance, bins=bins)
    return (
        annotated.groupby(["sample_id", "clonal_rank_bin"], observed=False, sort=False)
        .agg(
            n_clonotypes=("clone_id", "nunique"),
            occupied_frequency=("frequency", "sum"),
            count=("count", "sum"),
        )
        .reset_index()
    )


def public_repertoire(
    abundance: pd.DataFrame,
    *,
    min_samples: int = 2,
    min_count: float = 0,
    min_frequency: float = 0,
) -> pd.DataFrame:
    """Create a cohort-level table of clonotypes shared across samples."""

    if min_samples < 1:
        raise ValueError("min_samples must be positive")
    table = _with_frequency(abundance)
    table = table[(table["count"] >= min_count) & (table["frequency"] >= min_frequency)]
    metadata_columns = [
        column
        for column in ("locus", "junction", "junction_aa", "v_call", "d_call", "j_call", "c_call")
        if column in table
    ]
    aggregations: dict[str, tuple[str, object]] = {
        "n_samples": ("sample_id", "nunique"),
        "sample_ids": ("sample_id", lambda values: ";".join(sorted(set(map(str, values))))),
        "total_count": ("count", "sum"),
        "mean_frequency": ("frequency", "mean"),
        "max_frequency": ("frequency", "max"),
    }
    aggregations.update({column: (column, "first") for column in metadata_columns})
    result = table.groupby("clone_id", observed=True).agg(**aggregations).reset_index()
    result = result[result["n_samples"] >= min_samples]
    return result.sort_values(
        ["n_samples", "total_count"], ascending=[False, False]
    ).reset_index(drop=True)


def public_overlap(
    abundance: pd.DataFrame,
    *,
    metric: str = "jaccard",
    min_count: float = 0,
    min_frequency: float = 0,
) -> pd.DataFrame:
    """Calculate intersection, Jaccard or overlap-coefficient similarity."""

    if metric not in {"intersection", "jaccard", "overlap"}:
        raise ValueError("metric must be intersection, jaccard or overlap")
    table = _with_frequency(abundance)
    table = table[(table["count"] >= min_count) & (table["frequency"] >= min_frequency)]
    presence = abundance_matrix(table).gt(0).astype(int)
    intersection = presence @ presence.T
    if metric == "intersection":
        return intersection.astype(float)
    richness = np.diag(intersection.to_numpy()).astype(float)
    if metric == "jaccard":
        denominator = richness[:, None] + richness[None, :] - intersection.to_numpy()
    else:
        denominator = np.minimum(richness[:, None], richness[None, :])
    values = np.divide(
        intersection.to_numpy(dtype=float),
        denominator,
        out=np.zeros_like(denominator, dtype=float),
        where=denominator > 0,
    )
    return pd.DataFrame(values, index=presence.index, columns=presence.index)
