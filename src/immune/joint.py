"""Clone-aware RNA/VDJ summaries, bulk matching, and expression integration."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import entropy, fisher_exact

from ._objects import airr_data, cell_obs, record, sync_obs
from ._optional import require_dependency
from .clonotypes import CloneDefinition, add_clonotype_ids, clone_abundance
from .stats import adjust_pvalues


def clone_summary(
    data,
    *,
    sample_col="sample_id",
    clone_key="clone_id",
    airr_mod="airr",
    gex_mod="gex",
    key_added="clone_summary",
):
    """Sample-specific clone sizes; never substitute UMI counts for cell counts."""
    cells = cell_obs(data, airr_mod=airr_mod, gex_mod=gex_mod, clone_key=clone_key)
    valid = cells.dropna(subset=[sample_col, clone_key])
    valid = valid[~valid[clone_key].astype(str).isin(["None", "nan", "NaN", ""])]
    table = (
        valid.groupby([sample_col, clone_key], observed=True).size().rename("count").reset_index()
    )
    table = table.rename(columns={sample_col: "sample_id", clone_key: "clone_id"})
    totals = table.groupby("sample_id", observed=True)["count"].transform("sum")
    table["frequency"] = table["count"] / totals
    table["count_unit"] = "cell"
    record(
        data,
        key_added,
        table,
        params={
            "denominator": "clone_annotated_airr_cells",
            "sample_col": sample_col,
            "clone_key": clone_key,
        },
    )
    return table


def clonal_expansion(
    data,
    *,
    sample_col="sample_id",
    clone_key="clone_id",
    min_cells=2,
    airr_mod="airr",
    gex_mod="gex",
    key_added="clonal_expansion",
):
    """Annotate per-sample clone size and expanded/singleton status."""
    if min_cells < 2:
        raise ValueError("min_cells must be at least two")
    table = clone_summary(
        data, sample_col=sample_col, clone_key=clone_key, airr_mod=airr_mod, gex_mod=gex_mod
    )
    adata = airr_data(data, airr_mod)
    keys = pd.MultiIndex.from_frame(adata.obs[[sample_col, clone_key]])
    sizes = table.set_index(["sample_id", "clone_id"])["count"].reindex(keys)
    adata.obs["clone_size"] = sizes.to_numpy()
    status = pd.Series(pd.NA, index=adata.obs_names, dtype="string")
    detected = adata.obs["clone_size"].notna()
    status.loc[detected] = np.where(
        adata.obs.loc[detected, "clone_size"].ge(min_cells), "expanded", "not_expanded"
    )
    adata.obs[key_added] = status
    sync_obs(data, airr_mod)
    return data


def phenotype_diversity(
    data, *, composition_key="phenotype_composition", key_added="phenotype_diversity"
):
    from .get import result

    composition = data if isinstance(data, pd.DataFrame) else result(data, composition_key)
    rows = []
    for (sample, clone), group in composition.groupby(["sample_id", "clone_id"], observed=True):
        probabilities = group["within_clone_fraction"].to_numpy()
        shannon = float(entropy(probabilities))
        rows.append(
            {
                "sample_id": sample,
                "clone_id": clone,
                "n_states": int((probabilities > 0).sum()),
                "shannon": shannon,
                "effective_states": float(np.exp(shannon)),
                "n_cells": int(group["cell_count"].sum()),
            }
        )
    table = pd.DataFrame(
        rows,
        columns=["sample_id", "clone_id", "n_states", "shannon", "effective_states", "n_cells"],
    )
    if not isinstance(data, pd.DataFrame):
        record(data, key_added, table, params={"entropy_base": "e"}, backend="scipy")
    return table


def clone_state_enrichment(
    data,
    *,
    composition_key="phenotype_composition",
    min_cells=3,
    key_added="clone_state_enrichment",
):
    """Exploratory within-sample Fisher tests, not biological-replicate inference."""
    from .get import result

    composition = data if isinstance(data, pd.DataFrame) else result(data, composition_key)
    rows = []
    for sample, table in composition.groupby("sample_id", observed=True):
        matrix = table.pivot(index="clone_id", columns="phenotype", values="cell_count").fillna(0)
        total = float(matrix.to_numpy().sum())
        for clone, counts in matrix.iterrows():
            clone_total = float(counts.sum())
            if clone_total < min_cells:
                continue
            for phenotype, a in counts.items():
                b = clone_total - a
                c = matrix[phenotype].sum() - a
                d = total - a - b - c
                tested = fisher_exact([[int(a), int(b)], [int(c), int(d)]], alternative="greater")
                rows.append(
                    {
                        "sample_id": sample,
                        "clone_id": clone,
                        "phenotype": phenotype,
                        "odds_ratio": float(tested.statistic),
                        "p_value": float(tested.pvalue),
                        "method": "within_sample_fisher",
                    }
                )
    table = pd.DataFrame(
        rows, columns=["sample_id", "clone_id", "phenotype", "odds_ratio", "p_value", "method"]
    )
    table["fdr"] = adjust_pvalues(table["p_value"].tolist())
    if not isinstance(data, pd.DataFrame):
        record(
            data,
            key_added,
            table,
            params={"min_cells": min_cells, "fdr_family": "all_clone_state_tests"},
            backend="scipy",
        )
    return table


def match_bulk(
    data,
    bulk_chains,
    *,
    sample_pairs=None,
    definition=None,
    clone_key="clone_id",
    airr_mod="airr",
    key_added="bulk_matches",
):
    """Match a single locus to bulk while preserving full paired clone identities.

    The default compares equal sample IDs only. Cross-time/tissue comparisons
    require an explicit bulk_sample_id/sc_sample_id mapping. Donor metadata,
    when present on both inputs, is checked and never silently crossed.
    """
    from .get import chain_table

    definition = definition or CloneDefinition()
    if len(definition.loci) != 1:
        raise ValueError("Bulk matching requires a single-locus definition")
    adata = airr_data(data, airr_mod)
    if clone_key not in adata.obs:
        raise KeyError(clone_key)
    chains = chain_table(data, airr_mod=airr_mod)
    cells = add_clonotype_ids(chains, definition, key_added="match_key")
    if "productive" in cells:
        cells = cells[cells["productive"].fillna(False)]
    cells = cells.dropna(subset=["match_key"])
    cells["paired_clone_id"] = cells["cell_id"].map(adata.obs[clone_key])
    abundance = clone_abundance(bulk_chains, definition).rename(
        columns={
            "sample_id": "bulk_sample_id",
            "clone_id": "match_key",
            "count": "bulk_count",
            "frequency": "bulk_frequency",
            "count_unit": "bulk_count_unit",
        }
    )
    if sample_pairs is None:
        shared = sorted(set(abundance["bulk_sample_id"]) & set(cells["sample_id"]))
        sample_pairs = pd.DataFrame({"bulk_sample_id": shared, "sc_sample_id": shared})
    pairs = sample_pairs[["bulk_sample_id", "sc_sample_id"]].drop_duplicates()
    if pairs.isna().any().any():
        raise ValueError("Sample pairs must not contain missing IDs")
    if set(pairs["bulk_sample_id"]) - set(abundance["bulk_sample_id"]):
        raise ValueError("Unknown bulk sample in sample_pairs")
    if set(pairs["sc_sample_id"]) - set(adata.obs["sample_id"]):
        raise ValueError("Unknown single-cell sample in sample_pairs")
    if "donor_id" in bulk_chains and "donor_id" in adata.obs:
        for frame in (bulk_chains, adata.obs):
            if (
                frame.groupby("sample_id", observed=True)["donor_id"]
                .nunique(dropna=False)
                .gt(1)
                .any()
            ):
                raise ValueError("A sample must have exactly one donor")
        bd = bulk_chains.drop_duplicates("sample_id").set_index("sample_id")["donor_id"]
        sd = adata.obs.drop_duplicates("sample_id").set_index("sample_id")["donor_id"]
        for pair in pairs.itertuples(index=False):
            if pd.isna(bd[pair.bulk_sample_id]) or pd.isna(sd[pair.sc_sample_id]):
                raise ValueError("Donor identifiers must not be missing")
            if bd[pair.bulk_sample_id] != sd[pair.sc_sample_id]:
                raise ValueError("Bulk matching cannot cross donors")
    left = cells[["sample_id", "cell_id", "paired_clone_id", "match_key"]].drop_duplicates()
    left = left.rename(columns={"sample_id": "sc_sample_id"}).merge(pairs, on="sc_sample_id")
    table = left.merge(
        abundance[
            ["bulk_sample_id", "match_key", "bulk_count", "bulk_frequency", "bulk_count_unit"]
        ],
        on=["bulk_sample_id", "match_key"],
        how="left",
    )
    table["matched"] = table["bulk_count"].notna()
    table["n_paired_clones"] = table.groupby(["sc_sample_id", "match_key"], observed=True)[
        "paired_clone_id"
    ].transform("nunique")
    table["ambiguous"] = table["n_paired_clones"].gt(1)
    table["clone_definition"] = definition.label
    record(
        data,
        key_added,
        table,
        params={
            "definition": definition.label,
            "matching_locus": definition.loci[0],
            "sample_pairs": pairs.to_dict(orient="list"),
        },
    )
    return table


def annotate_bulk_matches(data, *, result_key="bulk_matches", airr_mod="airr", prefix="bulk"):
    """Annotate a single bulk comparison per cell; reject ambiguous annotations."""
    from .get import result

    table = result(data, result_key)
    if table["cell_id"].duplicated().any():
        raise ValueError(
            "Select one bulk comparison and one matching chain per cell before annotation"
        )
    adata = airr_data(data, airr_mod)
    indexed = table.set_index("cell_id")
    for col in ("matched", "ambiguous", "bulk_count", "bulk_frequency", "bulk_count_unit"):
        adata.obs[f"{prefix}_{col.removeprefix('bulk_')}"] = indexed[col].reindex(adata.obs_names)
    sync_obs(data, airr_mod)
    return data


def clone_pseudobulk(
    data,
    *,
    groupby="clonal_expansion",
    sample_col="sample_id",
    layer="counts",
    airr_mod="airr",
    gex_mod="gex",
    metadata_cols=(),
    **kwargs,
):
    """Delegate count aggregation to cellscope using receptor-defined groups."""
    cs = require_dependency("cellscope", extra="joint", feature="Clone-aware expression analysis")
    if not hasattr(data, "mod"):
        raise TypeError("Clone-aware pseudobulk requires MuData with gex and airr")
    airr = airr_data(data, airr_mod)
    gex = data.mod[gex_mod]
    labels = airr.obs[groupby].reindex(gex.obs_names)
    subset = gex[labels.notna()].copy()
    subset.obs[groupby] = labels.reindex(subset.obs_names)
    return cs.tl.pseudobulk(
        subset,
        sample_col=sample_col,
        groups_col=groupby,
        layer=layer,
        metadata_cols=metadata_cols,
        **kwargs,
    )


def clone_expression(
    data,
    *,
    design,
    contrast,
    method="pydeseq2",
    groupby="clonal_expansion",
    metadata_cols=(),
    min_cells=10,
    pseudobulk_kwargs=None,
    de_kwargs=None,
    key_added="clone_expression",
):
    """Compare conditions within receptor-defined groups using sample-level DE."""
    cs = require_dependency(
        "cellscope", extra="joint", feature="Clone-aware differential expression"
    )
    pdata = clone_pseudobulk(
        data,
        groupby=groupby,
        metadata_cols=metadata_cols,
        min_cells=min_cells,
        **(pseudobulk_kwargs or {}),
    )
    options = dict(de_kwargs or {})
    options.setdefault("method", method)
    options.setdefault("sample_col", (pseudobulk_kwargs or {}).get("sample_col", "sample_id"))
    table, models = cs.tl.differential_expression(
        pdata, design=design, contrast=contrast, groups_col=groupby, **options
    )
    selected_method = str(table["method"].iloc[0])
    record(
        data,
        key_added,
        table,
        params={
            "design": design,
            "contrast": list(contrast),
            "groupby": groupby,
            "replicate_unit": "sample",
            "method": selected_method,
        },
        backend="edgepython" if selected_method.startswith("edgepython") else selected_method,
    )
    return table, pdata, models
