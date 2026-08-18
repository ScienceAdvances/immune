"""Optional Matplotlib visualizations for immune."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from .stats import timecourse


def _plot_modules():
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError("Install immune[plot] to use plotting functions") from exc
    return plt, sns


def plot_clone_trajectories(
    abundance: pd.DataFrame,
    clones: Sequence[str],
    *,
    sample_order: Sequence[str] | None = None,
    ax: object | None = None,
    log_y: bool = True,
) -> object:
    """Plot frequency trajectories for selected clonotypes."""

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - depends on optional extra
        raise ImportError("Install immune[plot] to use plotting functions") from exc
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 4))
    matrix = timecourse(abundance, sample_order=sample_order, value="frequency")
    for clone_id in clones:
        if clone_id not in matrix.index:
            continue
        ax.plot(matrix.columns, matrix.loc[clone_id], marker="o", label=clone_id)
    if log_y:
        ax.set_yscale("log")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Clonotype frequency")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    return ax


def plot_bulk_diversity(
    diversity: pd.DataFrame,
    *,
    metric: str | None = None,
    hue: str | None = None,
    ax: object | None = None,
) -> object:
    """Plot a scikit-bio diversity result using seaborn."""

    plt, sns = _plot_modules()
    data = diversity.copy()
    if {"metric", "value"}.issubset(data.columns):
        if metric is not None:
            data = data[data["metric"] == metric]
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 4))
        sns.barplot(data=data, x="sample_id", y="value", hue=hue, ax=ax)
        ax.set_ylabel(metric or "Diversity")
    else:
        if metric is None or metric not in data.columns:
            raise KeyError("For wide summaries, metric must name a result column")
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 4))
        sns.barplot(data=data, x="sample_id", y=metric, hue=hue, ax=ax)
    ax.tick_params(axis="x", rotation=45)
    return ax


def plot_segment_usage(
    usage: pd.DataFrame,
    *,
    ax: object | None = None,
    cmap: str = "viridis",
    **kwargs: object,
) -> object:
    """Plot a sample-by-V/D/J matrix as a seaborn heatmap."""

    plt, sns = _plot_modules()
    if ax is None:
        width = max(7.0, usage.shape[1] * 0.28)
        _, ax = plt.subplots(figsize=(width, max(3.0, usage.shape[0] * 0.45)))
    sns.heatmap(usage, cmap=cmap, ax=ax, **kwargs)
    ax.set_xlabel(usage.columns.name or "Segment")
    ax.set_ylabel("Sample")
    return ax


def plot_spectratype(
    data: pd.DataFrame,
    *,
    value: str = "frequency",
    ax: object | None = None,
    **kwargs: object,
) -> object:
    """Plot Python-native CDR3 length spectra."""

    plt, sns = _plot_modules()
    if value not in data:
        raise KeyError(value)
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4))
    sns.lineplot(
        data=data,
        x="length",
        y=value,
        hue="sample_id",
        marker="o",
        ax=ax,
        **kwargs,
    )
    ax.set_xlabel("CDR3 length")
    return ax


def plot_repertoire_overlap(
    distances: pd.DataFrame,
    *,
    ax: object | None = None,
    cmap: str = "mako",
    **kwargs: object,
) -> object:
    """Plot a scikit-bio repertoire distance matrix."""

    plt, sns = _plot_modules()
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(distances, cmap=cmap, square=True, ax=ax, **kwargs)
    return ax


def plot_hill_diversity(
    profile: pd.DataFrame,
    *,
    ax: object | None = None,
    **kwargs: object,
) -> object:
    """Plot Hill number profiles across diversity orders."""

    plt, sns = _plot_modules()
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 4))
    sns.lineplot(
        data=profile,
        x="q",
        y="hill_number",
        hue="sample_id",
        marker="o",
        ax=ax,
        **kwargs,
    )
    ax.set_xlabel("Hill order (q)")
    ax.set_ylabel("Hill number")
    return ax


def plot_coverage_diversity(
    coverage: pd.DataFrame,
    *,
    relative: bool = False,
    ax: object | None = None,
    **kwargs: object,
) -> object:
    """Plot DXX coverage diversity for one or more percentages."""

    plt, sns = _plot_modules()
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 4))
    y = "dxx_fraction" if relative else "dxx"
    sns.lineplot(
        data=coverage,
        x="percentage",
        y=y,
        hue="sample_id",
        marker="o",
        ax=ax,
        **kwargs,
    )
    ax.set_ylabel("Fraction of clonotypes" if relative else "Number of clonotypes")
    return ax


def plot_clonality(
    clonality: pd.DataFrame,
    *,
    bin_col: str | None = None,
    value: str = "occupied_frequency",
    ax: object | None = None,
    **kwargs: object,
) -> object:
    """Plot occupied repertoire space as stacked clonality bins."""

    plt, _ = _plot_modules()
    if bin_col is None:
        candidates = [
            column
            for column in ("clonal_prop_bin", "clonal_rank_bin")
            if column in clonality
        ]
        if not candidates:
            raise KeyError("A clonal_prop_bin or clonal_rank_bin column is required")
        bin_col = candidates[0]
    if value not in clonality:
        raise KeyError(value)
    matrix = clonality.pivot_table(
        index="sample_id",
        columns=bin_col,
        values=value,
        aggfunc="sum",
        fill_value=0,
        observed=False,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4))
    matrix.plot(kind="bar", stacked=True, ax=ax, **kwargs)
    ax.set_xlabel("Sample")
    ax.set_ylabel(value.replace("_", " ").title())
    ax.legend(title=bin_col, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    return ax


def plot_rank_abundance(
    ranked: pd.DataFrame,
    *,
    log_x: bool = True,
    log_y: bool = True,
    ax: object | None = None,
    **kwargs: object,
) -> object:
    """Plot clone rank-abundance curves."""

    plt, sns = _plot_modules()
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 4))
    sns.lineplot(
        data=ranked,
        x="rank",
        y="frequency",
        hue="sample_id",
        ax=ax,
        **kwargs,
    )
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    return ax


def plot_public_repertoire(
    public: pd.DataFrame,
    *,
    top_n: int = 30,
    value: str = "n_samples",
    ax: object | None = None,
    **kwargs: object,
) -> object:
    """Plot the most widely shared public clonotypes."""

    plt, sns = _plot_modules()
    if value not in public:
        raise KeyError(value)
    data = public.nlargest(top_n, [value, "total_count"])
    if ax is None:
        _, ax = plt.subplots(figsize=(8, max(4, len(data) * 0.25)))
    sns.barplot(data=data, x=value, y="clone_id", ax=ax, **kwargs)
    return ax
