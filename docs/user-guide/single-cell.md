# Single-cell RNA and V(D)J

`immune` provides thin adapters around Scirpy and Scanpy rather than
reimplementing their receptor and transcriptome algorithms.

## Convert receptor chains to Scirpy

Start from cell-level Cell Ranger contigs:

```python
chains = iu.io.read_10x(
    "cellranger/outs/filtered_contig_annotations.csv",
    sample_id="sample_01",
)
airr = iu.pp.to_scirpy(chains)
```

All qualifying chains are retained. When the same barcode occurs in multiple
samples, cell IDs become `sample_id:barcode`; otherwise the original barcode is
preserved for direct transcriptome alignment.

## Merge with gene expression

```python
import scanpy as sc

adata = sc.read_h5ad("rna.h5ad")
mdata = iu.pp.merge_with_transcriptome(adata, airr)
```

The result is a MuData object with gene-expression and AIRR modalities. Verify
that cell identifiers and sample annotations agree before analysis.

## Chain QC and clonotypes

```python
iu.pp.scirpy_qc(mdata)

iu.tl.scirpy_define_clonotypes(
    mdata,
    receptor_arms="all",
    dual_ir="primary_only",
    same_v_gene=False,
    same_j_gene=False,
)
```

Arguments are passed to the installed Scirpy version. Set receptor-arm and
dual-receptor policies explicitly in reproducible studies.

## Repertoire summaries

```python
iu.tl.scirpy_summary(mdata, groupby="sample_id")

iu.pl.clonal_expansion(mdata, groupby="sample_id")
iu.pl.repertoire_overlap(mdata, groupby="sample_id")
iu.pl.vdj_usage(mdata, full_combination=False)
iu.pl.clonotype_network(mdata, color="cell_type")
```

For a convenient default orchestration:

```python
iu.tl.scirpy_repertoire(mdata, groupby="sample_id")
```

The decomposed QC, clonotyping and summary functions are preferred when a study
must document or change scientific parameters.

## Transcriptome preprocessing

```python
iu.tl.scanpy_workflow(
    adata,
    layer="counts",
    n_top_genes=3000,
    n_neighbors=15,
)
```

This calls a conventional Scanpy preprocessing/embedding sequence. For an
established project, use the study's existing Scanpy workflow and pass the
annotated object to the receptor integration steps.

## scVI and neighborhood differential abundance

```python
model = iu.tl.scvi(
    adata,
    layer="counts",
    batch_key="sample_id",
)

mdata_milo, milo = iu.tl.milo(
    adata,
    sample_col="sample_id",
    design="~ condition",
    neighbors_kwargs={"use_rep": "X_scVI"},
)
```

The native scvi-tools model and Pertpy Milo objects are returned. Keep their
version and fitted parameters in the study provenance.

## Pseudobulk expression

```python
pseudobulk = iu.tl.pseudobulk(
    adata,
    sample_col="sample_id",
    groups_col="cell_type",
    layer="counts",
)
```

Pseudobulk aggregation delegates to decoupler. Use sample-level replicates for
downstream expression inference rather than treating cells as independent
biological replicates.
