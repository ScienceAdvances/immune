"""Scirpy adapters for receptor specificity, BCR maturation and clone structure."""

from __future__ import annotations

from ._objects import airr_data, record, sync_obs
from ._optional import require_dependency


def receptor_query(
    data,
    reference,
    *,
    sequence="aa",
    metric="identity",
    cutoff=None,
    include_ref_cols=None,
    strategy="unique-only",
    same_v_gene=False,
    receptor_arms="all",
    dual_ir="any",
    airr_mod="airr",
    airr_mod_ref="airr",
    key_added="specificity",
    distance_kwargs=None,
    query_kwargs=None,
):
    """Match receptors to an explicit Scirpy reference and publish its annotations.

    Database similarity is a specificity hypothesis. Keep ambiguous matches
    explicit with strategy='json', or unique-only (default). Junction/CDR3
    conventions and donor HLA restrictions must be harmonized before querying.
    Returns the reference-annotation DataFrame in query-cell order.
    """
    ir = require_dependency("scirpy", extra="singlecell", feature="Receptor specificity query")
    adata, ref = airr_data(data, airr_mod), airr_data(reference, airr_mod_ref)
    ir.pp.index_chains(adata)
    ir.pp.index_chains(ref)
    distance_key, query_key = f"{key_added}_dist", f"{key_added}_query"
    options = dict(distance_kwargs or {})
    if cutoff is not None:
        options["cutoff"] = cutoff
    ir.pp.ir_dist(adata, ref, sequence=sequence, metric=metric, key_added=distance_key, **options)
    ir.tl.ir_query(
        adata,
        ref,
        sequence=sequence,
        metric=metric,
        distance_key=distance_key,
        key_added=query_key,
        same_v_gene=same_v_gene,
        receptor_arms=receptor_arms,
        dual_ir=dual_ir,
        **(query_kwargs or {}),
    )
    table = ir.tl.ir_query_annotate(
        adata,
        ref,
        sequence=sequence,
        metric=metric,
        query_key=query_key,
        include_ref_cols=include_ref_cols,
        strategy=strategy,
        inplace=False,
    )
    for col in table:
        adata.obs[f"{key_added}_{col}"] = table[col].reindex(adata.obs_names)
    record(
        data,
        key_added,
        table,
        params={
            "sequence": sequence,
            "metric": str(metric),
            "cutoff": cutoff,
            "strategy": strategy,
            "same_v_gene": same_v_gene,
            "receptor_arms": receptor_arms,
            "dual_ir": dual_ir,
        },
        backend="scirpy",
    )
    sync_obs(data, airr_mod)
    return table


def bcr_clonotypes(
    data,
    *,
    cutoff=15,
    airr_mod="airr",
    scope="donor_id",
    key_added="bcr_clone_id",
    distance_kwargs=None,
    clonotype_kwargs=None,
):
    """Cluster somatically mutated BCR junctions by normalized nucleotide Hamming distance."""
    from .singlecell import scirpy_define_clonotypes

    return scirpy_define_clonotypes(
        data,
        sequence="nt",
        metric="normalized_hamming",
        airr_mod=airr_mod,
        scope=scope,
        key_added=key_added,
        distance_kwargs={"cutoff": cutoff, **(distance_kwargs or {})},
        clonotype_kwargs=clonotype_kwargs,
    )


def mutational_load(data, *, airr_mod="airr", regions=("full",), **kwargs):
    """Compute BCR mutation burden from AIRR sequence/germline alignments.

    Region-specific estimates additionally require AIRR region coordinate
    fields. 10x junction-only CSVs do not contain the required germline alignments.
    """
    ir = require_dependency("scirpy", extra="singlecell", feature="BCR mutational load")
    adata = airr_data(data, airr_mod)
    fields = set(adata.obsm["airr"].fields)
    ak = require_dependency("awkward", extra="singlecell", feature="AIRR alignment validation")
    for key in (
        kwargs.get("sequence_key", "sequence_alignment"),
        kwargs.get("germline_key", "germline_alignment_d_mask"),
    ):
        if key not in fields or not any(ak.to_list(ak.flatten(adata.obsm["airr"][key], axis=1))):
            raise ValueError(
                f"Mutational load requires AIRR field {key!r}; import aligned AIRR data"
            )
    ir.pp.index_chains(adata)
    ir.tl.mutational_load(adata, regions=list(regions), **kwargs)
    record(data, "mutational_load", params={"regions": list(regions)}, backend="scirpy")
    sync_obs(data, airr_mod)
    return data


def clonotype_network(data, *, airr_mod="airr", **kwargs):
    """Build a Scirpy clonotype graph/layout for immune.pl.clonotype_network."""
    ir = require_dependency("scirpy", extra="singlecell", feature="Clonotype network")
    airr_data(data, airr_mod)
    ir.tl.clonotype_network(data, airr_mod=airr_mod, **kwargs)
    record(data, "clonotype_network", params=kwargs, backend="scirpy")
    return data


def clonotype_modularity(
    data, *, target_col="clone_id", connectivity_key="gex:connectivities", airr_mod="airr", **kwargs
):
    """Test whether clonotypes occupy localized transcriptomic graph neighborhoods."""
    ir = require_dependency("scirpy", extra="singlecell", feature="Clonotype modularity")
    airr_data(data, airr_mod)
    ir.tl.clonotype_modularity(
        data, target_col=target_col, connectivity_key=connectivity_key, airr_mod=airr_mod, **kwargs
    )
    record(
        data,
        "clonotype_modularity",
        params={"target_col": target_col, "connectivity_key": connectivity_key, **kwargs},
        backend="scirpy",
    )
    sync_obs(data, airr_mod)
    return data
