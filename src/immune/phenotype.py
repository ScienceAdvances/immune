"""Clone-aware single-cell phenotype summaries and PhenoTrack-style flow."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .clonotypes import CloneDefinition, add_clonotype_ids


def _cell_clone_table(
    chains: pd.DataFrame,
    definition: CloneDefinition,
) -> pd.DataFrame:
    assigned = add_clonotype_ids(chains, definition)
    assigned = assigned.dropna(subset=["cell_id", "clone_id"]).copy()
    if assigned.empty:
        return pd.DataFrame(columns=["sample_id", "cell_id", "clone_id"])
    score = (
        assigned["umi_count"]
        .fillna(assigned["read_count"])
        .fillna(0)
        .astype(float)
    )
    assigned["_support"] = score
    # Dual beta chains are retained in the chain table. For cell annotation,
    # select the best-supported chain and leave dual-chain modeling to a later
    # explicit analysis rather than silently counting a cell twice.
    assigned = assigned.sort_values("_support", ascending=False).drop_duplicates(
        ["sample_id", "cell_id"]
    )
    return assigned[["sample_id", "cell_id", "clone_id"]]


def phenotype_composition(
    single_cell_chains: pd.DataFrame,
    cell_metadata: pd.DataFrame,
    *,
    phenotype_col: str,
    definition: CloneDefinition | None = None,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
) -> pd.DataFrame:
    """Count phenotypes within each clonotype and sample.

    ``cell_metadata`` may be ``adata.obs`` or an ordinary DataFrame. If
    ``cell_id_col`` is absent, its index is interpreted as the cell identifier.
    """

    definition = definition or CloneDefinition()
    if phenotype_col not in cell_metadata:
        raise KeyError(phenotype_col)
    metadata = cell_metadata.copy()
    if cell_id_col not in metadata:
        metadata[cell_id_col] = metadata.index.astype(str)
    if sample_col not in metadata:
        samples = single_cell_chains["sample_id"].dropna().unique()
        if len(samples) != 1:
            raise KeyError(
                f"{sample_col!r} is required when receptor data contains multiple samples"
            )
        metadata[sample_col] = str(samples[0])

    cell_clones = _cell_clone_table(single_cell_chains, definition)
    merged = cell_clones.merge(
        metadata[[sample_col, cell_id_col, phenotype_col]],
        left_on=["sample_id", "cell_id"],
        right_on=[sample_col, cell_id_col],
        how="inner",
    )
    if merged.empty:
        return pd.DataFrame(
            columns=[
                "sample_id",
                "clone_id",
                "phenotype",
                "cell_count",
                "clone_cell_count",
                "sample_cell_count",
                "within_clone_fraction",
                "sample_fraction",
            ]
        )
    grouped = (
        merged.groupby(["sample_id", "clone_id", phenotype_col], observed=True)
        .size()
        .rename("cell_count")
        .reset_index()
        .rename(columns={phenotype_col: "phenotype"})
    )
    grouped["clone_cell_count"] = grouped.groupby(
        ["sample_id", "clone_id"], observed=True
    )["cell_count"].transform("sum")
    grouped["sample_cell_count"] = grouped.groupby("sample_id", observed=True)[
        "cell_count"
    ].transform("sum")
    grouped["within_clone_fraction"] = grouped["cell_count"] / grouped["clone_cell_count"]
    grouped["sample_fraction"] = grouped["cell_count"] / grouped["sample_cell_count"]
    return grouped.sort_values(
        ["sample_id", "clone_cell_count", "clone_id", "phenotype"],
        ascending=[True, False, True, True],
    ).reset_index(drop=True)


def phenotypic_flux(
    composition: pd.DataFrame,
    *,
    from_sample: str,
    to_sample: str,
    min_cells: int = 1,
) -> pd.DataFrame:
    """Compute the L1 change in within-clone phenotype composition."""

    matrix = composition.pivot_table(
        index=["sample_id", "clone_id"],
        columns="phenotype",
        values="cell_count",
        aggfunc="sum",
        fill_value=0,
    )
    if from_sample not in matrix.index.get_level_values("sample_id"):
        raise KeyError(from_sample)
    if to_sample not in matrix.index.get_level_values("sample_id"):
        raise KeyError(to_sample)
    start = matrix.xs(from_sample, level="sample_id")
    end = matrix.xs(to_sample, level="sample_id")
    clones = start.index.intersection(end.index)
    rows: list[dict[str, object]] = []
    for clone_id in clones:
        start_counts = start.loc[clone_id].to_numpy(dtype=float)
        end_counts = end.loc[clone_id].to_numpy(dtype=float)
        if start_counts.sum() < min_cells or end_counts.sum() < min_cells:
            continue
        start_prob = start_counts / start_counts.sum()
        end_prob = end_counts / end_counts.sum()
        rows.append(
            {
                "clone_id": clone_id,
                "from_sample": from_sample,
                "to_sample": to_sample,
                "from_cells": int(start_counts.sum()),
                "to_cells": int(end_counts.sum()),
                "phenotypic_flux": float(np.abs(start_prob - end_prob).sum()),
            }
        )
    result = pd.DataFrame.from_records(rows)
    if result.empty:
        return pd.DataFrame(
            columns=[
                "clone_id",
                "from_sample",
                "to_sample",
                "from_cells",
                "to_cells",
                "phenotypic_flux",
            ]
        )
    return result.sort_values("phenotypic_flux", ascending=False, ignore_index=True)


def phenotype_flow(
    composition: pd.DataFrame,
    *,
    from_sample: str,
    to_sample: str,
    clones: list[str] | None = None,
) -> pd.DataFrame:
    """Compute PhenoTrack's independent-approximation outflow and inflow.

    The result is a clone-level probabilistic flow between sampled cell-state
    distributions, not a lineage observation of the same physical cells.
    """

    phenotypes = sorted(composition["phenotype"].dropna().astype(str).unique())
    start = composition[composition["sample_id"].astype(str) == str(from_sample)]
    end = composition[composition["sample_id"].astype(str) == str(to_sample)]
    start_clones = set(start["clone_id"])
    end_clones = set(end["clone_id"])
    selected = start_clones & end_clones
    if clones is not None:
        selected &= set(clones)
    rows: list[dict[str, object]] = []
    for clone_id in sorted(selected):
        start_clone = start[start["clone_id"] == clone_id].set_index("phenotype")
        end_clone = end[end["clone_id"] == clone_id].set_index("phenotype")
        start_probs = np.array(
            [float(start_clone["within_clone_fraction"].get(p, 0)) for p in phenotypes]
        )
        end_probs = np.array(
            [float(end_clone["within_clone_fraction"].get(p, 0)) for p in phenotypes]
        )
        start_frequency = float(start_clone["sample_fraction"].sum())
        end_frequency = float(end_clone["sample_fraction"].sum())
        for i, from_phenotype in enumerate(phenotypes):
            for j, to_phenotype in enumerate(phenotypes):
                base = start_probs[i] * end_probs[j]
                rows.append(
                    {
                        "clone_id": clone_id,
                        "from_sample": from_sample,
                        "to_sample": to_sample,
                        "from_phenotype": from_phenotype,
                        "to_phenotype": to_phenotype,
                        "outflow": base * start_frequency,
                        "inflow": base * end_frequency,
                    }
                )
    return pd.DataFrame.from_records(rows)
