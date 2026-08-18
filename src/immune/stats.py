"""Repertoire summaries and longitudinal expansion statistics."""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
import pandas as pd

try:  # SciPy is a package dependency; fallback keeps lightweight source tests usable.
    from scipy.stats import fisher_exact as _scipy_fisher_exact
except ImportError:  # pragma: no cover - exercised only in minimal development runtimes
    _scipy_fisher_exact = None


def _hypergeom_probability(x: int, row1: int, row2: int, col1: int) -> float:
    return math.comb(row1, x) * math.comb(row2, col1 - x) / math.comb(row1 + row2, col1)


def fisher_exact_two_sided(a: float, b: float, c: float, d: float) -> float:
    """Two-sided Fisher exact test with a small-table dependency-free fallback."""

    cells = [round(max(0.0, float(value))) for value in (a, b, c, d)]
    a_i, b_i, c_i, d_i = cells
    if _scipy_fisher_exact is not None:
        result = _scipy_fisher_exact([[a_i, b_i], [c_i, d_i]], alternative="two-sided")
        return float(getattr(result, "pvalue", result[1]))
    row1, row2 = a_i + b_i, c_i + d_i
    col1 = a_i + c_i
    total = row1 + row2
    if total == 0:
        return 1.0
    if total > 100_000:
        raise ImportError("SciPy is required for Fisher tests with more than 100,000 counts")
    low = max(0, col1 - row2)
    high = min(row1, col1)
    observed = _hypergeom_probability(a_i, row1, row2, col1)
    p_value = sum(
        _hypergeom_probability(x, row1, row2, col1)
        for x in range(low, high + 1)
        if _hypergeom_probability(x, row1, row2, col1) <= observed + 1e-12
    )
    return min(1.0, float(p_value))


def adjust_pvalues(pvalues: Sequence[float], method: str = "bh") -> np.ndarray:
    """Adjust p-values using Benjamini-Hochberg or Bonferroni."""

    values = np.asarray(pvalues, dtype=float)
    if values.size == 0:
        return values
    if method.lower() in {"bonferroni", "bonf"}:
        return np.clip(values * values.size, 0, 1)
    if method.lower() not in {"bh", "fdr_bh", "benjamini-hochberg"}:
        raise ValueError("method must be 'bh' or 'bonferroni'")
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * values.size / np.arange(1, values.size + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.clip(adjusted, 0, 1)
    return output


def timecourse(
    abundance: pd.DataFrame,
    *,
    sample_order: Sequence[str] | None = None,
    value: str = "frequency",
    pseudocount: float | None = None,
) -> pd.DataFrame:
    """Return a clonotype-by-sample trajectory matrix."""

    if value not in abundance:
        raise KeyError(value)
    matrix = abundance.pivot_table(
        index="clone_id", columns="sample_id", values=value, aggfunc="sum", fill_value=0
    )
    if sample_order is not None:
        matrix = matrix.reindex(columns=list(sample_order), fill_value=0)
    if pseudocount is not None:
        matrix = matrix.mask(matrix == 0, pseudocount)
    return matrix


def _expansion_pvalue(
    baseline_count: float,
    baseline_total: float,
    comparison_count: float,
    comparison_total: float,
    fold_change: float,
) -> float:
    if baseline_total <= 0 or comparison_total <= 0:
        return 1.0
    baseline_frequency = baseline_count / baseline_total
    comparison_frequency = comparison_count / comparison_total
    if comparison_frequency < fold_change * baseline_frequency:
        return 1.0
    effective_baseline_total = max(round(baseline_total / fold_change), 1)
    baseline_int = round(baseline_count)
    comparison_int = round(comparison_count)
    comparison_total_int = max(round(comparison_total), comparison_int)
    if baseline_int > effective_baseline_total:
        return 1.0
    return fisher_exact_two_sided(
        baseline_int,
        effective_baseline_total - baseline_int,
        comparison_int,
        comparison_total_int - comparison_int,
    )


def test_longitudinal_expansion(
    abundance: pd.DataFrame,
    *,
    baseline_sample: str,
    comparison_samples: Sequence[str] | None = None,
    fold_change: float = 2.0,
    correction: str = "bh",
    pseudocount: float = 1 / 3,
    de_novo_count_max: float | None = None,
) -> pd.DataFrame:
    """Identify clones expanding after a baseline sample.

    This keeps CloneTrack's modified Fisher concept but reports effect sizes,
    the best time point and FDR. MiXCR fractional clone counts are rounded only
    for the exact test; original counts and frequencies remain in the output.
    """

    required = {"sample_id", "clone_id", "count"}
    missing = required - set(abundance.columns)
    if missing:
        raise KeyError(f"abundance is missing columns: {sorted(missing)}")
    samples = list(dict.fromkeys(abundance["sample_id"].astype(str)))
    if baseline_sample not in samples:
        raise KeyError(f"Unknown baseline sample: {baseline_sample}")
    if comparison_samples is None:
        baseline_index = samples.index(baseline_sample)
        comparison_samples = samples[baseline_index + 1 :]
    comparison_samples = list(comparison_samples)
    if not comparison_samples:
        raise ValueError("At least one comparison sample is required")
    unknown = set(comparison_samples) - set(samples)
    if unknown:
        raise KeyError(f"Unknown comparison samples: {sorted(unknown)}")

    count_matrix = abundance.pivot_table(
        index="clone_id", columns="sample_id", values="count", aggfunc="sum", fill_value=0
    )
    count_matrix = count_matrix.reindex(
        columns=[baseline_sample, *comparison_samples], fill_value=0
    )
    totals = count_matrix.sum(axis=0)
    rows: list[dict[str, object]] = []
    for clone_id, counts in count_matrix.iterrows():
        baseline_count = float(counts[baseline_sample])
        if de_novo_count_max is not None and baseline_count > de_novo_count_max:
            continue
        baseline_frequency = (baseline_count + pseudocount) / (
            float(totals[baseline_sample]) + pseudocount
        )
        candidates: list[dict[str, object]] = []
        for sample in comparison_samples:
            comparison_count = float(counts[sample])
            comparison_frequency = (comparison_count + pseudocount) / (
                float(totals[sample]) + pseudocount
            )
            candidates.append(
                {
                    "sample": sample,
                    "comparison_count": comparison_count,
                    "comparison_frequency": comparison_frequency,
                    "fold_change": comparison_frequency / baseline_frequency,
                    "p_value": _expansion_pvalue(
                        baseline_count,
                        float(totals[baseline_sample]),
                        comparison_count,
                        float(totals[sample]),
                        fold_change,
                    ),
                }
            )
        best = min(candidates, key=lambda item: (float(item["p_value"]), -float(item["fold_change"])))
        time_adjusted_p = min(1.0, float(best["p_value"]) * len(comparison_samples))
        rows.append(
            {
                "clone_id": clone_id,
                "baseline_sample": baseline_sample,
                "best_sample": best["sample"],
                "baseline_count": baseline_count,
                "comparison_count": best["comparison_count"],
                "baseline_frequency": baseline_frequency,
                "comparison_frequency": best["comparison_frequency"],
                "fold_change": best["fold_change"],
                "p_value": best["p_value"],
                "time_adjusted_p_value": time_adjusted_p,
            }
        )
    result = pd.DataFrame.from_records(rows)
    if result.empty:
        return result
    result["q_value"] = adjust_pvalues(result["time_adjusted_p_value"], correction)
    result["expanded"] = (result["fold_change"] >= fold_change) & (result["q_value"] < 0.05)
    return result.sort_values(["q_value", "fold_change"], ascending=[True, False]).reset_index(
        drop=True
    )


# This is a public statistical test, not a pytest test case.  The marker keeps
# pytest from collecting it when users import the function into a test module.
test_longitudinal_expansion.__test__ = False


def repertoire_metrics(abundance: pd.DataFrame) -> pd.DataFrame:
    """Compute richness, Shannon entropy, Simpson diversity and clonality."""

    records: list[dict[str, object]] = []
    for sample_id, sample in abundance.groupby("sample_id", observed=True):
        counts = sample["count"].to_numpy(dtype=float)
        counts = counts[counts > 0]
        total = counts.sum()
        probabilities = counts / total if total > 0 else np.array([], dtype=float)
        richness = len(probabilities)
        shannon = float(-(probabilities * np.log(probabilities)).sum()) if richness else 0.0
        simpson = float(1 - np.square(probabilities).sum()) if richness else 0.0
        clonality = float(1 - shannon / np.log(richness)) if richness > 1 else 1.0 if richness else 0.0
        records.append(
            {
                "sample_id": sample_id,
                "total_count": total,
                "richness": richness,
                "shannon": shannon,
                "simpson": simpson,
                "clonality": clonality,
            }
        )
    return pd.DataFrame.from_records(records)
