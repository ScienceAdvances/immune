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
    score = assigned["umi_count"].fillna(assigned["read_count"]).fillna(0).astype(float)
    assigned["_support"] = score
    # Dual beta chains are retained in the chain table. For cell annotation,
    # select the best-supported chain and leave dual-chain modeling to a later
    # explicit analysis rather than silently counting a cell twice.
    assigned = assigned.sort_values("_support", ascending=False).drop_duplicates(
        ["sample_id", "cell_id"]
    )
    return assigned[["sample_id", "cell_id", "clone_id"]]


def _phenotype_composition_table(
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
    grouped["clone_cell_count"] = grouped.groupby(["sample_id", "clone_id"], observed=True)[
        "cell_count"
    ].transform("sum")
    grouped["sample_cell_count"] = grouped.groupby("sample_id", observed=True)[
        "cell_count"
    ].transform("sum")
    grouped["within_clone_fraction"] = grouped["cell_count"] / grouped["clone_cell_count"]
    grouped["sample_fraction"] = grouped["cell_count"] / grouped["sample_cell_count"]
    return grouped.sort_values(
        ["sample_id", "clone_cell_count", "clone_id", "phenotype"],
        ascending=[True, False, True, True],
    ).reset_index(drop=True)


def _phenotypic_flux_table(
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


def _phenotype_flow_table(
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


def phenotype_composition(
    data,
    cell_metadata=None,
    *,
    phenotype_col="cell_state",
    definition=None,
    cell_id_col="cell_id",
    sample_col="sample_id",
    clone_key="clone_id",
    airr_mod="airr",
    gex_mod="gex",
    key_added="phenotype_composition",
):
    """Count clone/state composition from native objects or legacy chain tables.

    For native objects, clonotypes come from the existing receptor annotation.
    The denominator is cells with both a clone and the requested state label.
    Original receptor records and RNA expression matrices are not modified.
    """
    if isinstance(data, pd.DataFrame):
        if cell_metadata is None:
            raise ValueError("cell_metadata is required for chain-table input")
        return _phenotype_composition_table(
            data,
            cell_metadata,
            phenotype_col=phenotype_col,
            definition=definition,
            cell_id_col=cell_id_col,
            sample_col=sample_col,
        )
    from ._objects import cell_obs, record

    cells = cell_obs(data, airr_mod=airr_mod, gex_mod=gex_mod, clone_key=clone_key)
    if definition is not None:
        raise ValueError("Define native clonotypes first; definition only applies to table input")
    columns = [sample_col, clone_key, phenotype_col]
    valid = cells.dropna(subset=columns)
    valid = valid[~valid[clone_key].astype(str).isin(["None", "nan", "NaN", ""])]
    result = valid.groupby(columns, observed=True).size().rename("cell_count").reset_index()
    result = result.rename(
        columns={sample_col: "sample_id", clone_key: "clone_id", phenotype_col: "phenotype"}
    )
    result["clone_cell_count"] = result.groupby(["sample_id", "clone_id"], observed=True)[
        "cell_count"
    ].transform("sum")
    result["sample_cell_count"] = result.groupby("sample_id", observed=True)[
        "cell_count"
    ].transform("sum")
    result["within_clone_fraction"] = result["cell_count"] / result["clone_cell_count"]
    result["sample_fraction"] = result["cell_count"] / result["sample_cell_count"]
    record(
        data,
        key_added,
        result,
        params={
            "phenotype_col": phenotype_col,
            "clone_key": clone_key,
            "denominator": "clone_and_state_annotated_cells",
            "n_input_cells": len(cells),
            "n_used_cells": len(valid),
        },
    )
    return result


def phenotypic_flux(
    data,
    *,
    from_sample,
    to_sample,
    min_cells=1,
    composition_key="phenotype_composition",
    key_added="phenotypic_flux",
):
    from ._objects import record
    from .get import result

    composition = data if isinstance(data, pd.DataFrame) else result(data, composition_key)
    table = _phenotypic_flux_table(
        composition, from_sample=from_sample, to_sample=to_sample, min_cells=min_cells
    )
    if not isinstance(data, pd.DataFrame):
        record(
            data,
            key_added,
            table,
            params={
                "from_sample": from_sample,
                "to_sample": to_sample,
                "min_cells": min_cells,
                "metric": "L1",
            },
        )
    return table


def phenotype_flow(
    data,
    *,
    from_sample,
    to_sample,
    clones=None,
    composition_key="phenotype_composition",
    key_added="phenotype_flow",
):
    """Independent-approximation flow of clone state distributions, not lineage."""
    from ._objects import record
    from .get import result

    composition = data if isinstance(data, pd.DataFrame) else result(data, composition_key)
    table = _phenotype_flow_table(
        composition, from_sample=from_sample, to_sample=to_sample, clones=clones
    )
    if not isinstance(data, pd.DataFrame):
        record(
            data,
            key_added,
            table,
            params={
                "from_sample": from_sample,
                "to_sample": to_sample,
                "model": "independent_state_distributions",
            },
        )
    return table
