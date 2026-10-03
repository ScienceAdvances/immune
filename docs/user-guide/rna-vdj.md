# RNA and VDJ joint analysis

Use cellscope for expression and cell-state analysis, and immune for
receptors, clones, phenotype distributions, and bulk repertoire links.
Both operate on the same native MuData object.

```python
import cellscope as cs
import immune as iu

mdata = iu.datasets.toy_multimodal()
cs.pp.qc(mdata, min_genes=5)
cs.tl.workflow(mdata, n_top_genes=30, n_pcs=10)

iu.tl.define_clonotypes(mdata, scope="donor_id")
iu.tl.clonal_expansion(mdata)
composition = iu.tl.phenotype_composition(mdata, phenotype_col="cell_state")
flux = iu.tl.phenotypic_flux(mdata, from_sample="s0", to_sample="s1")
flow = iu.tl.phenotype_flow(mdata, from_sample="s0", to_sample="s1")
tracking = iu.tl.track_clones(mdata)
expansion = iu.tl.longitudinal_expansion(mdata)

links = iu.tl.match_bulk(mdata, iu.datasets.toy_bulk())
iu.tl.annotate_bulk_matches(mdata)
iu.pl.clone_embedding(mdata, color="airr:clonal_expansion")
iu.pl.phenotype_flow(mdata)
iu.io.write(mdata, "study.h5mu")
```

For real paired 10x libraries, use the same library identifier on both sides:

```python
gex = cs.io.read_10x("filtered_feature_bc_matrix.h5", library_id="lib1",
                     sample_id="post", donor_id="patient1")
airr = iu.io.read_10x_vdj("all_contig_annotations.json", library_id="lib1",
                         sample_id="post", donor_id="patient1")
mdata = iu.pp.merge_with_transcriptome(gex, airr)
```

Bulk links use single-locus keys while full paired clonotypes remain
separate. Inspect `ambiguous` and `n_paired_clones` before interpreting a
TRB match. Supply explicit `sample_pairs` when linking different timepoints.
Phenotype flow is a distribution estimate, not cell lineage observation.

For sample-level expression inference, use raw counts and a biological
design. The toy dataset demonstrates mechanics; actual studies require
adequate donors and sampling.

```python
rna = mdata.mod["gex"]
labels = mdata.mod["airr"].obs["clonal_expansion"].reindex(rna.obs_names)
selected = rna[labels.notna()].copy()
selected.obs["clonal_expansion"] = labels.reindex(selected.obs_names)
pdata = cs.tl.pseudobulk(
    selected, groups_col="clonal_expansion", min_cells=3,
    metadata_cols=["condition", "donor_id"]
)
de, models = cs.tl.differential_expression(
    pdata, groups_col="clonal_expansion", design="~ donor_id + condition",
    contrast=["condition", "post", "pre"]
)
```
