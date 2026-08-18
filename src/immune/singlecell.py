"""Scirpy adapters for canonical single-cell immune receptor data."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import pandas as pd

from ._optional import require_dependency
from .schema import validate_rearrangements


def _python_value(value: object) -> object:
    if value is None or pd.isna(value):
        return None
    return value.item() if hasattr(value, "item") else value


def to_scirpy(
    chains: pd.DataFrame,
    *,
    make_cell_ids_unique: bool = True,
    separator: str = ":",
):
    """Convert the canonical chain table to a Scirpy AIRR AnnData object."""

    ir = require_dependency("scirpy", extra="singlecell", feature="Scirpy conversion")
    validate_rearrangements(chains)
    cell_chains = chains[chains["cell_id"].notna()].copy()
    if cell_chains.empty:
        raise ValueError("Scirpy conversion requires cell-level receptor chains")

    duplicated_barcodes = cell_chains.groupby("cell_id", observed=True)["sample_id"].nunique().gt(1)
    duplicated = set(duplicated_barcodes[duplicated_barcodes].index.astype(str))
    airr_cells = []
    for (sample_id, cell_id), cell_table in cell_chains.groupby(
        ["sample_id", "cell_id"], observed=True, sort=False
    ):
        original_cell_id = str(cell_id)
        use_composite = make_cell_ids_unique and original_cell_id in duplicated
        scirpy_cell_id = (
            f"{sample_id}{separator}{original_cell_id}" if use_composite else original_cell_id
        )
        cell = ir.io.AirrCell(scirpy_cell_id)
        cell["sample_id"] = str(sample_id)
        cell["original_cell_id"] = original_cell_id
        for _, row in cell_table.iterrows():
            chain = {
                "sequence_id": _python_value(row.get("sequence_id")),
                "locus": _python_value(row.get("locus")),
                "productive": _python_value(row.get("productive")),
                "junction": _python_value(row.get("junction")),
                "junction_aa": _python_value(row.get("junction_aa")),
                "v_call": _python_value(row.get("v_call")),
                "d_call": _python_value(row.get("d_call")),
                "j_call": _python_value(row.get("j_call")),
                "c_call": _python_value(row.get("c_call")),
                "consensus_count": _python_value(row.get("read_count")),
                "duplicate_count": _python_value(row.get("umi_count")),
            }
            cell.add_chain(chain)
        airr_cells.append(cell)
    return ir.io.from_airr_cells(airr_cells)


def merge_with_transcriptome(rna, airr, *, rna_mod: str = "gex", airr_mod: str = "airr"):
    """Merge gene-expression and AIRR AnnData objects in a MuData container."""

    mudata = require_dependency("mudata", extra="singlecell", feature="RNA/AIRR integration")
    return mudata.MuData({rna_mod: rna, airr_mod: airr})


def scirpy_qc(data, **kwargs: Any):
    """Index chains and run Scirpy's receptor quality-control annotations."""

    ir = require_dependency("scirpy", extra="singlecell", feature="Scirpy quality control")
    ir.pp.index_chains(data, **kwargs.pop("index_chains", {}))
    ir.tl.chain_qc(data, **kwargs)
    return data


def scirpy_define_clonotypes(
    data,
    *,
    sequence: str = "nt",
    metric: str = "identity",
    distance_kwargs: dict[str, Any] | None = None,
    clonotype_kwargs: dict[str, Any] | None = None,
):
    """Compute receptor distances and define clonotypes with Scirpy."""

    ir = require_dependency("scirpy", extra="singlecell", feature="Scirpy clonotypes")
    ir.pp.ir_dist(data, sequence=sequence, metric=metric, **(distance_kwargs or {}))
    ir.tl.define_clonotypes(data, **(clonotype_kwargs or {}))
    return data


def scirpy_summary(
    data,
    *,
    groupby: str,
    target_col: str = "clone_id",
    diversity_metrics: Iterable[str] = ("normalized_shannon_entropy", "D50"),
    overlap_measure: str = "jaccard",
):
    """Add Scirpy clonal expansion, diversity, overlap and spectratype summaries."""

    ir = require_dependency("scirpy", extra="singlecell", feature="Scirpy summaries")
    ir.tl.clonal_expansion(data, target_col=target_col, expanded_in=groupby)
    for metric in diversity_metrics:
        ir.tl.alpha_diversity(
            data,
            groupby,
            target_col=target_col,
            metric=metric,
            key_added=f"alpha_diversity_{metric}",
        )
    ir.tl.repertoire_overlap(
        data,
        groupby,
        target_col=target_col,
        overlap_measure=overlap_measure,
    )
    spectratype = ir.tl.spectratype(data, target_col=groupby)
    if hasattr(data, "uns"):
        data.uns[f"spectratype_{groupby}"] = spectratype
    return data


def run_scirpy_repertoire(
    data,
    *,
    groupby: str,
    sequence: str = "nt",
    metric: str = "identity",
    distance_kwargs: dict[str, Any] | None = None,
    clonotype_kwargs: dict[str, Any] | None = None,
):
    """Run the recommended Scirpy QC-to-summary repertoire workflow."""

    scirpy_qc(data)
    scirpy_define_clonotypes(
        data,
        sequence=sequence,
        metric=metric,
        distance_kwargs=distance_kwargs,
        clonotype_kwargs=clonotype_kwargs,
    )
    return scirpy_summary(data, groupby=groupby)


def plot_scirpy(data, kind: str, **kwargs: Any):
    """Call a named, supported Scirpy plotting function."""

    ir = require_dependency("scirpy", extra="singlecell", feature="Scirpy plotting")
    supported = {
        "alpha_diversity",
        "clonal_expansion",
        "clonotype_imbalance",
        "clonotype_modularity",
        "clonotype_network",
        "embedding",
        "group_abundance",
        "logoplot_cdr3_motif",
        "repertoire_overlap",
        "spectratype",
        "vdj_usage",
    }
    if kind not in supported:
        raise ValueError(f"Unsupported Scirpy plot {kind!r}; choose from {sorted(supported)}")
    return getattr(ir.pl, kind)(data, **kwargs)
