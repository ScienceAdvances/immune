"""Donor-separated detection trajectories and CloneTrack-style expansion."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ._objects import record
from .joint import clone_summary
from .stats import adjust_pvalues, test_longitudinal_expansion


def _inputs(data, sample_metadata, donor_col, time_col):
    abundance = data.copy() if isinstance(data, pd.DataFrame) else clone_summary(data)
    if sample_metadata is None:
        if isinstance(data, pd.DataFrame):
            raise ValueError("Sample metadata is required for abundance-table input")
        sample_metadata = data.obs[["sample_id", donor_col, time_col]].dropna().drop_duplicates()
    metadata = sample_metadata.copy()
    if "sample_id" not in metadata:
        metadata = metadata.rename_axis("sample_id").reset_index()
    if metadata[["sample_id", donor_col, time_col]].isna().any().any():
        raise ValueError("Sample, donor, and time metadata must not be missing")
    if metadata["sample_id"].duplicated().any():
        raise ValueError("Sample metadata must have one row per sample")
    if metadata.duplicated([donor_col, time_col]).any():
        raise ValueError("Use one sample per donor and timepoint")
    time = metadata[time_col]
    ordered_category = isinstance(time.dtype, pd.CategoricalDtype) and time.dtype.ordered
    if not (
        pd.api.types.is_numeric_dtype(time)
        or pd.api.types.is_datetime64_any_dtype(time)
        or ordered_category
    ):
        raise ValueError("Time must be numeric, datetime, or an ordered categorical")
    if set(abundance["sample_id"]) - set(metadata["sample_id"]):
        raise ValueError("Some abundance samples have no study metadata")
    if not np.isfinite(abundance["count"]).all() or abundance["count"].lt(0).any():
        raise ValueError("Abundance counts must be finite and nonnegative")
    return abundance, metadata.sort_values([donor_col, time_col])


def track_clones(
    data,
    sample_metadata=None,
    *,
    donor_col="donor_id",
    time_col="timepoint",
    min_count=1,
    key_added="clone_tracking",
):
    """Complete observed clone trajectories within donors, including zero counts.

    A zero denotes not detected at the specified threshold, not biological
    extinction. Time values must be numeric, datetimes, or an ordered category.
    """
    if min_count <= 0:
        raise ValueError("min_count must be positive")
    abundance, metadata = _inputs(data, sample_metadata, donor_col, time_col)
    rows = []
    for donor, samples in metadata.groupby(donor_col, observed=True, sort=False):
        sample_ids = samples["sample_id"].tolist()
        selected = abundance[abundance["sample_id"].isin(sample_ids)]
        counts = selected.pivot_table(
            index="clone_id", columns="sample_id", values="count", aggfunc="sum", fill_value=0
        ).reindex(columns=sample_ids, fill_value=0)
        totals = counts.sum(axis=0)
        for clone_id, series in counts.iterrows():
            detected = series.ge(min_count)
            for sample_id, count in series.items():
                rows.append(
                    {
                        "donor_id": donor,
                        "sample_id": sample_id,
                        "clone_id": clone_id,
                        "timepoint": samples.set_index("sample_id").loc[sample_id, time_col],
                        "count": count,
                        "frequency": count / totals[sample_id] if totals[sample_id] else np.nan,
                        "detected": bool(detected[sample_id]),
                        "n_detected_timepoints": int(detected.sum()),
                        "persistent": bool(detected.all()),
                    }
                )
    table = pd.DataFrame(
        rows,
        columns=[
            "donor_id",
            "sample_id",
            "clone_id",
            "timepoint",
            "count",
            "frequency",
            "detected",
            "n_detected_timepoints",
            "persistent",
        ],
    )
    if not isinstance(data, pd.DataFrame):
        record(
            data,
            key_added,
            table,
            params={"min_count": min_count, "time_col": time_col, "donor_col": donor_col},
        )
    return table


def longitudinal_expansion(
    data,
    *,
    baseline_sample=None,
    comparison_samples=None,
    sample_metadata=None,
    donor_col="donor_id",
    time_col="timepoint",
    key_added="longitudinal_expansion",
    **kwargs,
):
    """Legacy table tests or donor-specific native-object longitudinal tests.

    Native mode uses the earliest timepoint per donor, or a donor->sample
    baseline mapping. Bonferroni corrects timepoint selection; BH correction
    then covers all tested donor/clone combinations.
    """
    if isinstance(data, pd.DataFrame) and sample_metadata is None:
        if baseline_sample is None:
            raise ValueError("baseline_sample is required for legacy table input")
        return test_longitudinal_expansion(
            data, baseline_sample=baseline_sample, comparison_samples=comparison_samples, **kwargs
        )
    abundance, metadata = _inputs(data, sample_metadata, donor_col, time_col)
    if kwargs.get("correction", "bh") != "bh":
        raise ValueError("Native study mode uses global BH correction")
    if comparison_samples is not None:
        raise ValueError("Select study samples in sample_metadata for donor-aware mode")
    if isinstance(baseline_sample, str) and metadata[donor_col].nunique() != 1:
        raise ValueError("Supply a donor->baseline mapping for multi-donor studies")
    results = []
    for donor, samples in metadata.groupby(donor_col, observed=True, sort=False):
        ids = samples["sample_id"].tolist()
        baseline = (
            baseline_sample.get(donor) if isinstance(baseline_sample, dict) else baseline_sample
        )
        baseline = baseline or ids[0]
        if baseline not in ids:
            raise ValueError(f"Baseline sample does not belong to donor {donor!r}")
        comparisons = ids[ids.index(baseline) + 1 :]
        if not comparisons:
            raise ValueError(f"No follow-up sample after baseline for donor {donor!r}")
        selected = abundance[abundance["sample_id"].isin([baseline, *comparisons])]
        if (
            selected.groupby("sample_id")["count"]
            .sum()
            .reindex([baseline, *comparisons], fill_value=0)
            .le(0)
            .any()
        ):
            raise ValueError(
                "Expansion testing requires nonzero sampled counts at every tested timepoint"
            )
        if "count_unit" in selected and selected["count_unit"].nunique() > 1:
            raise ValueError("Cannot compare different count units")
        table = test_longitudinal_expansion(
            selected, baseline_sample=baseline, comparison_samples=comparisons, **kwargs
        )
        if not table.empty:
            table["donor_id"] = donor
            results.append(table)
    result = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
    if not result.empty:
        result["q_value"] = adjust_pvalues(result["time_adjusted_p_value"].tolist())
        result["expanded"] = result["fold_change"].ge(kwargs.get("fold_change", 2.0)) & result[
            "q_value"
        ].lt(0.05)
    if not isinstance(data, pd.DataFrame):
        record(
            data,
            key_added,
            result,
            params={
                "donor_col": donor_col,
                "time_col": time_col,
                "fdr_family": "all_donor_clone_tests",
                **kwargs,
            },
        )
    return result
