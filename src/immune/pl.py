"""Plotting: Python-native bulk plots and Scirpy plotting adapters."""

from __future__ import annotations

from typing import Any

import pandas as pd

from .plotting import (
    plot_bulk_diversity,
    plot_clonality,
    plot_clone_trajectories,
    plot_coverage_diversity,
    plot_hill_diversity,
    plot_public_repertoire,
    plot_rank_abundance,
    plot_repertoire_overlap,
    plot_segment_usage,
    plot_spectratype,
)
from .singlecell import plot_scirpy

# Scanpy-style plotting names omit the ``plot_`` prefix inside ``iu.pl``.
bulk_diversity = plot_bulk_diversity
clone_trajectories = plot_clone_trajectories
segment_usage = plot_segment_usage
scirpy = plot_scirpy
bulk_repertoire_overlap = plot_repertoire_overlap
bulk_spectratype = plot_spectratype
clonality = plot_clonality
coverage_diversity = plot_coverage_diversity
hill_diversity = plot_hill_diversity
public_overlap = plot_repertoire_overlap
public_repertoire = plot_public_repertoire
rank_abundance = plot_rank_abundance


def alpha_diversity(data, **kwargs: Any):
    return plot_scirpy(data, "alpha_diversity", **kwargs)


def clonal_expansion(data, **kwargs: Any):
    return plot_scirpy(data, "clonal_expansion", **kwargs)


def clonotype_imbalance(data, **kwargs: Any):
    return plot_scirpy(data, "clonotype_imbalance", **kwargs)


def clonotype_modularity(data, **kwargs: Any):
    return plot_scirpy(data, "clonotype_modularity", **kwargs)


def clonotype_network(data, **kwargs: Any):
    return plot_scirpy(data, "clonotype_network", **kwargs)


def embedding(data, **kwargs: Any):
    return plot_scirpy(data, "embedding", **kwargs)


def group_abundance(data, **kwargs: Any):
    return plot_scirpy(data, "group_abundance", **kwargs)


def repertoire_overlap(data, **kwargs: Any):
    """Plot a bulk distance matrix or Scirpy repertoire overlap."""

    if isinstance(data, pd.DataFrame):
        return plot_repertoire_overlap(data, **kwargs)
    return plot_scirpy(data, "repertoire_overlap", **kwargs)


def spectratype(data, **kwargs: Any):
    """Plot a bulk length table or a Scirpy spectratype."""

    if isinstance(data, pd.DataFrame):
        return plot_spectratype(data, **kwargs)
    return plot_scirpy(data, "spectratype", **kwargs)


repertoire_overlap_sc = repertoire_overlap
spectratype_sc = spectratype


def vdj_usage(data, **kwargs: Any):
    return plot_scirpy(data, "vdj_usage", **kwargs)


__all__ = [
    "alpha_diversity",
    "bulk_diversity",
    "bulk_repertoire_overlap",
    "bulk_spectratype",
    "clonal_expansion",
    "clonality",
    "clone_trajectories",
    "clonotype_imbalance",
    "clonotype_modularity",
    "clonotype_network",
    "coverage_diversity",
    "embedding",
    "group_abundance",
    "hill_diversity",
    "plot_bulk_diversity",
    "plot_clonality",
    "plot_clone_trajectories",
    "plot_coverage_diversity",
    "plot_hill_diversity",
    "plot_public_repertoire",
    "plot_rank_abundance",
    "plot_repertoire_overlap",
    "plot_scirpy",
    "plot_segment_usage",
    "plot_spectratype",
    "public_overlap",
    "public_repertoire",
    "rank_abundance",
    "repertoire_overlap",
    "repertoire_overlap_sc",
    "scirpy",
    "segment_usage",
    "spectratype",
    "spectratype_sc",
    "vdj_usage",
]
