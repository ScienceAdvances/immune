"""Python-native bulk repertoire analysis backed by mature libraries."""

from __future__ import annotations

import re
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from ._optional import require_dependency

DEFAULT_ALPHA_METRICS = (
    "observed_features",
    "shannon",
    "simpson",
    "inv_simpson",
    "pielou_e",
    "gini_index",
    "chao1",
)


def abundance_matrix(
    abundance: pd.DataFrame,
    *,
    value_col: str = "count",
) -> pd.DataFrame:
    """Return a sample-by-clonotype matrix for mature diversity backends."""

    required = {"sample_id", "clone_id", value_col}
    missing = required.difference(abundance.columns)
    if missing:
        raise KeyError(f"Abundance table is missing columns: {sorted(missing)}")
    matrix = abundance.pivot_table(
        index="sample_id",
        columns="clone_id",
        values=value_col,
        aggfunc="sum",
        fill_value=0,
        observed=True,
    )
    matrix = matrix.astype(float)
    if (matrix.to_numpy() < 0).any():
        raise ValueError("Abundance values must be non-negative")
    return matrix


def skbio_alpha_diversity(
    abundance: pd.DataFrame,
    *,
    metrics: Sequence[str] = DEFAULT_ALPHA_METRICS,
    value_col: str = "count",
    round_counts: bool = True,
    **kwargs: Any,
) -> pd.DataFrame:
    """Compute alpha diversity through scikit-bio's validated metric drivers."""

    diversity = require_dependency(
        "skbio.diversity", extra="bulk", feature="scikit-bio alpha diversity"
    )
    matrix = abundance_matrix(abundance, value_col=value_col)
    if round_counts:
        matrix = matrix.round()
    pieces: list[pd.DataFrame] = []
    for metric in metrics:
        values = diversity.alpha_diversity(
            metric,
            matrix.to_numpy(),
            ids=matrix.index.astype(str).tolist(),
            **kwargs,
        )
        pieces.append(
            pd.DataFrame(
                {
                    "sample_id": values.index.astype(str),
                    "metric": metric,
                    "value": values.to_numpy(dtype=float),
                }
            )
        )
    if not pieces:
        return pd.DataFrame(columns=["sample_id", "metric", "value"])
    return pd.concat(pieces, ignore_index=True)


def skbio_beta_diversity(
    abundance: pd.DataFrame,
    *,
    metric: str = "braycurtis",
    value_col: str = "count",
    **kwargs: Any,
) -> pd.DataFrame:
    """Compute a sample distance matrix through scikit-bio."""

    diversity = require_dependency(
        "skbio.diversity", extra="bulk", feature="scikit-bio beta diversity"
    )
    matrix = abundance_matrix(abundance, value_col=value_col)
    distances = diversity.beta_diversity(
        metric,
        matrix.to_numpy(),
        ids=matrix.index.astype(str).tolist(),
        **kwargs,
    )
    return distances.to_data_frame()


def bulk_summary(
    abundance: pd.DataFrame,
    *,
    metrics: Sequence[str] = DEFAULT_ALPHA_METRICS,
    round_counts: bool = True,
) -> pd.DataFrame:
    """Build an immunarch-style sample summary using scikit-bio metrics."""

    required = {"sample_id", "clone_id", "count", "frequency"}
    missing = required.difference(abundance.columns)
    if missing:
        raise KeyError(f"Abundance table is missing columns: {sorted(missing)}")
    basic = (
        abundance.groupby("sample_id", observed=True)
        .agg(
            total_count=("count", "sum"),
            n_clonotypes=("clone_id", "nunique"),
            max_frequency=("frequency", "max"),
        )
        .reset_index()
    )
    diversity = skbio_alpha_diversity(
        abundance,
        metrics=metrics,
        round_counts=round_counts,
    ).pivot(index="sample_id", columns="metric", values="value")
    diversity.columns.name = None
    return basic.merge(diversity.reset_index(), on="sample_id", how="left")


def _gene_level(value: object, level: str, ambiguous: str) -> object:
    if value is None or pd.isna(value):
        return pd.NA
    calls = [item.strip() for item in re.split(r"[,;]", str(value)) if item.strip()]
    if not calls:
        return pd.NA
    if len(calls) > 1:
        if ambiguous == "exclude":
            return pd.NA
        if ambiguous == "first":
            calls = calls[:1]

    def transform(call: str) -> str:
        call = call.upper()
        if level in {"segment", "family"}:
            call = call.split("*", maxsplit=1)[0]
        if level == "family":
            call = call.split("-", maxsplit=1)[0]
        return call

    values = list(dict.fromkeys(transform(call) for call in calls))
    return ";".join(values)


def gene_usage(
    abundance: pd.DataFrame,
    *,
    gene: str = "v_call",
    level: str = "segment",
    ambiguous: str = "first",
    weighted: bool = True,
    normalize: bool = True,
    by_locus: bool = False,
    format: str = "wide",
) -> pd.DataFrame:
    """Summarize V/D/J/C usage by allele, segment or gene family."""

    if gene not in {"v_call", "d_call", "j_call", "c_call"}:
        raise ValueError("gene must be v_call, d_call, j_call or c_call")
    if level not in {"allele", "segment", "family"}:
        raise ValueError("level must be allele, segment or family")
    if ambiguous not in {"first", "exclude", "combine"}:
        raise ValueError("ambiguous must be first, exclude or combine")
    if format not in {"wide", "long"}:
        raise ValueError("format must be wide or long")
    required = {"sample_id", gene, "count"}
    if by_locus:
        required.add("locus")
    missing = required.difference(abundance.columns)
    if missing:
        raise KeyError(f"Abundance table is missing columns: {sorted(missing)}")
    table = abundance.copy()
    table["gene"] = table[gene].map(lambda value: _gene_level(value, level, ambiguous))
    table = table.dropna(subset=["gene"])
    table["_weight"] = table["count"] if weighted else 1.0
    groups = ["sample_id", *(["locus"] if by_locus else []), "gene"]
    result = table.groupby(groups, observed=True)["_weight"].sum().rename("usage").reset_index()
    if normalize:
        denominator_groups = ["sample_id", *(["locus"] if by_locus else [])]
        totals = result.groupby(denominator_groups, observed=True)["usage"].transform("sum")
        result["usage"] = np.divide(
            result["usage"],
            totals,
            out=np.zeros(len(result), dtype=float),
            where=totals.to_numpy() > 0,
        )
    result["gene_call"] = gene
    result["gene_level"] = level
    if format == "long":
        return result
    columns: str | list[str] = ["locus", "gene"] if by_locus else "gene"
    return result.pivot_table(
        index="sample_id",
        columns=columns,
        values="usage",
        aggfunc="sum",
        fill_value=0,
        observed=True,
    ).astype(float)


def segment_usage(
    abundance: pd.DataFrame,
    *,
    segment: str = "v_call",
    weighted: bool = True,
    normalize: bool = True,
    level: str = "segment",
    ambiguous: str = "first",
    by_locus: bool = False,
    format: str = "wide",
) -> pd.DataFrame:
    """Backward-compatible alias for :func:`gene_usage`."""

    return gene_usage(
        abundance,
        gene=segment,
        level=level,
        ambiguous=ambiguous,
        weighted=weighted,
        normalize=normalize,
        by_locus=by_locus,
        format=format,
    )


def spectratype(
    abundance: pd.DataFrame,
    *,
    sequence_col: str = "junction_aa",
    weighted: bool = True,
    normalize: bool = True,
) -> pd.DataFrame:
    """Summarize CDR3 length distributions in a tidy table."""

    if sequence_col not in {"junction", "junction_aa"}:
        raise ValueError("sequence_col must be junction or junction_aa")
    required = {"sample_id", sequence_col, "count"}
    missing = required.difference(abundance.columns)
    if missing:
        raise KeyError(f"Abundance table is missing columns: {sorted(missing)}")
    table = abundance.dropna(subset=[sequence_col]).copy()
    table["length"] = table[sequence_col].astype(str).str.len()
    table["abundance"] = table["count"] if weighted else 1.0
    result = (
        table.groupby(["sample_id", "length"], observed=True)["abundance"]
        .sum()
        .reset_index()
    )
    if normalize:
        totals = result.groupby("sample_id", observed=True)["abundance"].transform("sum")
        result["frequency"] = np.divide(
            result["abundance"],
            totals,
            out=np.zeros(len(result), dtype=float),
            where=totals.to_numpy() > 0,
        )
    return result


def rarefaction_curve(
    abundance: pd.DataFrame,
    *,
    depths: Sequence[int] | None = None,
    metric: str = "observed_features",
    n_iter: int = 20,
    seed: int | None = 0,
) -> pd.DataFrame:
    """Estimate rarefaction curves via scikit-bio count subsampling."""

    if n_iter < 1:
        raise ValueError("n_iter must be at least 1")
    diversity = require_dependency(
        "skbio.diversity", extra="bulk", feature="rarefaction diversity"
    )
    stats = require_dependency("skbio.stats", extra="bulk", feature="rarefaction subsampling")
    matrix = abundance_matrix(abundance).round().astype(int)
    rng = np.random.default_rng(seed)
    records: list[dict[str, object]] = []
    for sample_id, row in matrix.iterrows():
        counts = row.to_numpy(dtype=int)
        total = int(counts.sum())
        sample_depths = list(depths) if depths is not None else _default_depths(total)
        for depth in sample_depths:
            if depth < 1 or depth > total:
                continue
            estimates = []
            for _ in range(n_iter):
                subsampled = stats.subsample_counts(
                    counts,
                    int(depth),
                    replace=False,
                    seed=int(rng.integers(0, np.iinfo(np.int32).max)),
                )
                value = diversity.alpha_diversity(metric, subsampled)[0]
                estimates.append(float(value))
            values = np.asarray(estimates)
            records.append(
                {
                    "sample_id": str(sample_id),
                    "depth": int(depth),
                    "metric": metric,
                    "mean": float(values.mean()),
                    "std": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
                    "lower": float(np.quantile(values, 0.025)),
                    "upper": float(np.quantile(values, 0.975)),
                }
            )
    return pd.DataFrame.from_records(records)


def _default_depths(total: int, n_points: int = 20) -> list[int]:
    if total <= 0:
        return []
    return sorted(set(np.linspace(1, total, min(n_points, total), dtype=int).tolist()))


def differential_clonotype_abundance(
    abundance: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    *,
    design: str,
    contrast: Sequence[str],
    min_total_count: int = 10,
    min_samples: int = 2,
    alpha: float = 0.05,
    n_cpus: int | None = None,
    quiet: bool = False,
) -> tuple[pd.DataFrame, object, object]:
    """Fit replicate-aware clonotype differential abundance with PyDESeq2.

    Returns ``(results, dataset, statistics)`` so native PyDESeq2 plots and
    diagnostics remain accessible.
    """

    dds_module = require_dependency(
        "pydeseq2.dds", extra="differential", feature="clonotype differential abundance"
    )
    ds_module = require_dependency(
        "pydeseq2.ds", extra="differential", feature="clonotype differential abundance"
    )
    counts = abundance_matrix(abundance).round().astype(int)
    keep = (counts.sum(axis=0) >= min_total_count) & ((counts > 0).sum(axis=0) >= min_samples)
    counts = counts.loc[:, keep]
    if counts.empty:
        raise ValueError("No clonotypes pass the differential-abundance filters")
    if "sample_id" in sample_metadata.columns:
        metadata = sample_metadata.set_index("sample_id", drop=True).copy()
    else:
        metadata = sample_metadata.copy()
    metadata.index = metadata.index.astype(str)
    missing = counts.index.difference(metadata.index)
    if len(missing):
        raise KeyError(f"sample_metadata is missing samples: {missing.tolist()}")
    metadata = metadata.loc[counts.index]
    dds = dds_module.DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design=design,
        n_cpus=n_cpus,
        quiet=quiet,
    )
    dds.deseq2()
    statistics = ds_module.DeseqStats(
        dds,
        contrast=list(contrast),
        alpha=alpha,
        n_cpus=n_cpus,
        quiet=quiet,
    )
    statistics.summary()
    results = statistics.results_df.copy()
    results.index.name = "clone_id"
    return results.reset_index(), dds, statistics
