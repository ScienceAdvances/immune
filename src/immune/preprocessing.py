"""General repertoire filtering and count-preserving downsampling."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from ._optional import require_dependency


def filter_repertoire(
    table: pd.DataFrame,
    *,
    samples: Sequence[str] | None = None,
    loci: Sequence[str] | None = None,
    productive: bool | None = None,
    min_count: float | None = None,
    min_frequency: float | None = None,
    sequence_col: str | None = None,
    min_length: int | None = None,
    max_length: int | None = None,
    clone_ids: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Filter canonical chain or abundance tables with explicit criteria."""

    result = table.copy()
    if samples is not None:
        if "sample_id" not in result:
            raise KeyError("sample_id")
        result = result[result["sample_id"].astype(str).isin(map(str, samples))]
    if loci is not None:
        if "locus" not in result:
            raise KeyError("locus")
        wanted = {str(value).upper() for value in loci}
        result = result[result["locus"].astype(str).str.upper().isin(wanted)]
    if productive is not None:
        if "productive" not in result:
            raise KeyError("productive")
        result = result[result["productive"].astype("boolean").eq(productive)]
    if min_count is not None:
        if "count" not in result:
            raise KeyError("count")
        result = result[pd.to_numeric(result["count"], errors="coerce") >= min_count]
    if min_frequency is not None:
        if "frequency" not in result:
            raise KeyError("frequency")
        result = result[
            pd.to_numeric(result["frequency"], errors="coerce") >= min_frequency
        ]
    if clone_ids is not None:
        if "clone_id" not in result:
            raise KeyError("clone_id")
        result = result[result["clone_id"].astype(str).isin(map(str, clone_ids))]
    if min_length is not None or max_length is not None:
        sequence_col = sequence_col or "junction_aa"
        if sequence_col not in result:
            raise KeyError(sequence_col)
        lengths = result[sequence_col].astype("string").str.len()
        if min_length is not None:
            result = result[lengths >= min_length]
            lengths = lengths.loc[result.index]
        if max_length is not None:
            result = result[lengths <= max_length]
    return result.reset_index(drop=True)


def downsample_repertoire(
    abundance: pd.DataFrame,
    *,
    depth: int | str = "min",
    random_state: int | None = 0,
    drop_zero: bool = True,
) -> pd.DataFrame:
    """Downsample each repertoire without replacement using scikit-bio."""

    required = {"sample_id", "clone_id", "count"}
    missing = required.difference(abundance.columns)
    if missing:
        raise KeyError(f"Abundance table is missing columns: {sorted(missing)}")
    if abundance.duplicated(["sample_id", "clone_id"]).any():
        raise ValueError("Abundance must have one row per sample_id and clone_id")
    stats = require_dependency(
        "skbio.stats", extra="bulk", feature="repertoire downsampling"
    )
    table = abundance.copy()
    table["count"] = pd.to_numeric(table["count"], errors="raise").round().astype(int)
    totals = table.groupby("sample_id", observed=True)["count"].sum()
    if depth == "min":
        target = int(totals.min())
    elif isinstance(depth, int) and depth > 0:
        target = depth
    else:
        raise ValueError("depth must be a positive integer or 'min'")
    too_small = totals[totals < target]
    if not too_small.empty:
        raise ValueError(f"Samples below depth {target}: {too_small.index.tolist()}")

    rng = np.random.default_rng(random_state)
    pieces: list[pd.DataFrame] = []
    for _, sample in table.groupby("sample_id", observed=True, sort=False):
        sample = sample.copy()
        sample["original_count"] = sample["count"]
        sample["count"] = stats.subsample_counts(
            sample["count"].to_numpy(),
            target,
            replace=False,
            seed=int(rng.integers(0, np.iinfo(np.int32).max)),
        )
        if drop_zero:
            sample = sample[sample["count"] > 0]
        sample["frequency"] = sample["count"] / target if target else 0.0
        sample["downsample_depth"] = target
        pieces.append(sample)
    return pd.concat(pieces, ignore_index=True)
