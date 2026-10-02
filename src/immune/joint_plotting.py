"""Matplotlib joint-analysis plots; computations are explicit tools."""

from __future__ import annotations

import numpy as np
import pandas as pd


def _table(data, key):
    from .get import result

    return data if isinstance(data, pd.DataFrame) else result(data, key)


def phenotype_composition(
    data,
    *,
    key="phenotype_composition",
    sample_id=None,
    value="within_clone_fraction",
    ax=None,
    **kwargs,
):
    import matplotlib.pyplot as plt

    table = _table(data, key)
    if sample_id is not None:
        table = table[table["sample_id"].eq(sample_id)]
    matrix = table.pivot(index=["sample_id", "clone_id"], columns="phenotype", values=value).fillna(
        0
    )
    if ax is None:
        _, ax = plt.subplots()
    matrix.plot.bar(stacked=True, ax=ax, **kwargs)
    ax.set_ylabel(value)
    return ax


def phenotype_flow(data, *, key="phenotype_flow", value="outflow", ax=None, **kwargs):
    """Heatmap of the independent-approximation state flow."""
    import matplotlib.pyplot as plt

    table = _table(data, key)
    if table.empty:
        raise ValueError("No shared clones available for flow plotting")
    matrix = table.pivot_table(
        index="from_phenotype", columns="to_phenotype", values=value, aggfunc="sum", fill_value=0
    )
    if ax is None:
        _, ax = plt.subplots()
    artist = ax.imshow(matrix.to_numpy(), aspect="auto", **kwargs)
    ax.set_xticks(np.arange(len(matrix.columns)), matrix.columns, rotation=45)
    ax.set_yticks(np.arange(len(matrix.index)), matrix.index)
    ax.set(xlabel="State at follow-up", ylabel="State at baseline")
    ax.figure.colorbar(artist, ax=ax, label=value)
    return ax


def clone_embedding(
    data, *, color="airr:clonal_expansion", gex_mod="gex", basis="umap", ax=None, **kwargs
):
    """Overlay aligned receptor annotations on the RNA embedding."""
    from ._optional import require_dependency

    sc = require_dependency("scanpy", extra="singlecell", feature="Joint embedding plot")
    if not hasattr(data, "mod"):
        kwargs.setdefault("show", False)
        return sc.pl.embedding(data, basis=basis, color=color, ax=ax, **kwargs)
    gex = data.mod[gex_mod]
    if color not in gex.obs and color not in data.obs:
        raise KeyError(color)
    if color in gex.obs:
        kwargs.setdefault("show", False)
        return sc.pl.embedding(gex, basis=basis, color=color, ax=ax, **kwargs)
    # A temporary annotation is restored even if the backend raises.
    temporary = "_immune_plot_color"
    previous = gex.obs[temporary].copy() if temporary in gex.obs else None
    try:
        gex.obs[temporary] = data.obs[color].reindex(gex.obs_names)
        kwargs.setdefault("show", False)
        return sc.pl.embedding(gex, basis=basis, color=temporary, ax=ax, **kwargs)
    finally:
        if previous is None:
            del gex.obs[temporary]
        else:
            gex.obs[temporary] = previous
