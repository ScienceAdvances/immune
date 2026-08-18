"""Small adapters for AnnData/Scanpy without making them hard dependencies."""

from __future__ import annotations

import pandas as pd

from .clonotypes import CloneDefinition
from .phenotype import _cell_clone_table


def add_clonotypes_to_anndata(
    adata: object,
    single_cell_chains: pd.DataFrame,
    *,
    definition: CloneDefinition | None = None,
    sample_col: str | None = None,
    key_added: str = "immune_clone_id",
) -> None:
    """Add immune clonotypes to ``adata.obs`` in place.

    Cell barcodes are matched to ``adata.obs_names``. For multi-sample objects,
    pass the sample column in ``adata.obs`` so identical 10x barcodes from
    different libraries remain distinct.
    """

    if not hasattr(adata, "obs") or not hasattr(adata, "obs_names"):
        raise TypeError("adata must provide .obs and .obs_names like an AnnData object")
    definition = definition or CloneDefinition()
    cell_clones = _cell_clone_table(single_cell_chains, definition)
    obs = adata.obs
    obs_cells = pd.DataFrame({"_obs_name": adata.obs_names.astype(str)})
    obs_cells["cell_id"] = obs_cells["_obs_name"]
    if sample_col is None:
        samples = cell_clones["sample_id"].dropna().astype(str).unique()
        if len(samples) != 1:
            raise ValueError("sample_col is required for multi-sample receptor data")
        obs_cells["sample_id"] = samples[0]
    else:
        if sample_col not in obs:
            raise KeyError(sample_col)
        obs_cells["sample_id"] = obs[sample_col].astype(str).to_numpy()
    merged = obs_cells.merge(cell_clones, on=["sample_id", "cell_id"], how="left")
    obs[key_added] = pd.Series(
        merged["clone_id"].astype("string").to_numpy(), index=adata.obs_names, dtype="string"
    )
