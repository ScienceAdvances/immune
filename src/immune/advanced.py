"""Adapters to mature single-cell transcriptome and differential-abundance tools."""

from __future__ import annotations

from typing import Any

from ._optional import require_dependency


def scanpy_standard_workflow(
    adata,
    *,
    layer: str | None = None,
    n_top_genes: int = 2_000,
    n_neighbors: int = 15,
    n_pcs: int = 50,
    resolution: float = 1.0,
    copy: bool = False,
):
    """Run a conventional Scanpy normalization-to-UMAP workflow."""

    sc = require_dependency("scanpy", extra="singlecell", feature="Scanpy workflow")
    data = adata.copy() if copy else adata
    if layer is not None:
        data.X = data.layers[layer].copy()
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, n_top_genes=n_top_genes)
    sc.pp.pca(data, mask_var="highly_variable")
    sc.pp.neighbors(data, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(data)
    sc.tl.leiden(data, resolution=resolution)
    return data


def run_scvi(
    adata,
    *,
    layer: str | None = "counts",
    batch_key: str | None = None,
    categorical_covariate_keys: list[str] | None = None,
    model_kwargs: dict[str, Any] | None = None,
    train_kwargs: dict[str, Any] | None = None,
    key_added: str = "X_scVI",
):
    """Train scvi-tools SCVI and store its latent representation."""

    scvi = require_dependency("scvi", extra="advanced", feature="scvi-tools integration")
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
        categorical_covariate_keys=categorical_covariate_keys,
    )
    model = scvi.model.SCVI(adata, **(model_kwargs or {}))
    model.train(**(train_kwargs or {}))
    adata.obsm[key_added] = model.get_latent_representation()
    return model


def run_milo(
    adata,
    *,
    sample_col: str,
    design: str,
    prop: float = 0.1,
    neighbors_kwargs: dict[str, Any] | None = None,
    da_kwargs: dict[str, Any] | None = None,
):
    """Run Pertpy's Milo neighborhood differential-abundance workflow."""

    pt = require_dependency("pertpy", extra="advanced", feature="Milo differential abundance")
    sc = require_dependency("scanpy", extra="singlecell", feature="Milo neighbor graph")
    milo = pt.tl.Milo()
    mdata = milo.load(adata)
    sc.pp.neighbors(mdata["rna"], **(neighbors_kwargs or {}))
    milo.make_nhoods(mdata["rna"], prop=prop)
    mdata = milo.count_nhoods(mdata, sample_col=sample_col)
    milo.da_nhoods(mdata, design=design, **(da_kwargs or {}))
    milo.build_nhood_graph(mdata)
    return mdata, milo


def pseudobulk(
    adata,
    *,
    sample_col: str,
    groups_col: str,
    layer: str | None = None,
    mode: str = "sum",
    **kwargs: Any,
):
    """Aggregate cells using decoupler's pseudobulk implementation."""

    dc = require_dependency("decoupler", extra="advanced", feature="pseudobulk analysis")
    return dc.pp.pseudobulk(
        adata,
        sample_col=sample_col,
        groups_col=groups_col,
        layer=layer,
        mode=mode,
        **kwargs,
    )
