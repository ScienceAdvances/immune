# Single-cell V(D)J

`immune` provides receptor adapters around Scirpy. RNA algorithms belong to
cellscope; shared annotations and embeddings support receptor interpretation.

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
preserved for direct transcriptome alignment. Existing qualified `cell_id`
values remain intact when an original `barcode` column is present. An explicit
library identifier can qualify raw IDs. Final ID collisions raise an error;
the converter never invents numeric suffixes that would break RNA alignment.

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
    scope="donor_id",
    clonotype_kwargs={
        "receptor_arms": "all",
        "dual_ir": "primary_only",
        "same_v_gene": False,
        "same_j_gene": False,
    },
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

The combined workflow uses the output key from the selected sequence/metric
or `clonotype_kwargs["key_added"]` for all subsequent summaries, including
amino-acid runs. Choose `scope="donor_id"` explicitly for donor-restricted
clonotypes; its compatibility default `scope=None` groups shared sequences
across donors.

## Existing RNA annotations

Prepare the RNA object with cellscope in application code, then merge it with
AIRR data using `pp.merge_with_transcriptome`. See the RNA/VDJ guide. Receptor
plots and modularity analysis may use an existing RNA embedding or graph;
immune does not compute RNA preprocessing, embeddings, models or expression tests.
