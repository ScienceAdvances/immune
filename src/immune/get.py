"""Extract AIRR records, clone annotations, and stored joint-analysis tables."""

from __future__ import annotations

import copy

import pandas as pd

from ._objects import airr_data, cell_obs
from ._optional import require_dependency
from .schema import canonicalize


def obs_df(data, keys=None, *, airr_mod="airr", gex_mod="gex", clone_key="clone_id"):
    table = cell_obs(data, airr_mod=airr_mod, gex_mod=gex_mod, clone_key=clone_key)
    return (table if keys is None else table.loc[:, list(keys)]).copy()


def chain_table(data, *, airr_mod="airr"):
    """Flatten every stored chain, preserving AIRR fields and cell metadata."""
    ak = require_dependency("awkward", extra="singlecell", feature="AIRR extraction")
    adata = airr_data(data, airr_mod)
    rows = []
    for cell_id, chains in zip(adata.obs_names, ak.to_list(adata.obsm["airr"]), strict=True):
        metadata = adata.obs.loc[cell_id].to_dict()
        for chain in chains or []:
            row = {**chain, **metadata, "cell_id": str(cell_id)}
            row["sample_id"] = metadata.get("sample_id")
            row["read_count"] = chain.get("consensus_count")
            row["umi_count"] = chain.get("duplicate_count")
            row["source"] = "singlecell_airr"
            rows.append(row)
    return canonicalize(pd.DataFrame(rows))


def airr(data, fields, *, chain="VDJ_1", airr_mod="airr", **kwargs):
    """Retrieve selected indexed chains using Scirpy's public accessor."""
    ir = require_dependency("scirpy", extra="singlecell", feature="AIRR access")
    return ir.get.airr(data, fields, chain=chain, airr_mod=airr_mod, **kwargs)


def result(data, key):
    """Return a copy of a result table (or the full entry if no table is stored)."""
    entry = data.uns["immune"][key]
    return copy.deepcopy(entry.get("table", entry))
