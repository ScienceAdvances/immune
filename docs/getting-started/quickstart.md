# Quickstart

This workflow reads two bulk MiXCR repertoires and one Cell Ranger V(D)J
directory, then builds a shared clonotype definition for longitudinal and
single-cell integration.

## 1. Read receptor outputs

```python
import pandas as pd
import immune as iu

pre = iu.io.read_mixcr("data/pre.clones.tsv", sample_id="pre")
post = iu.io.read_mixcr("data/post.clones.tsv", sample_id="post")
cell_chains = iu.io.read_10x("data/cellranger/outs", sample_id="post_sc")

bulk_chains = pd.concat([pre, post], ignore_index=True)
```

Readers return the same AIRR-like chain schema regardless of source. Inspect
the parsed data before defining clonotypes:

```python
bulk_chains[[
    "sample_id", "locus", "junction", "junction_aa",
    "v_call", "j_call", "read_count", "umi_count",
]].head()
```

## 2. Define clonotypes explicitly

The conservative cross-platform default below uses nucleotide CDR3 and TRB,
while ignoring V/J calls that may differ between pipelines.

```python
definition = iu.pp.CloneDefinition(
    sequence="junction",
    loci=("TRB",),
)

bulk = iu.pp.clone_abundance(bulk_chains, definition)
```

Within one platform, a stricter definition can include V/J genes:

```python
strict = iu.pp.CloneDefinition(
    sequence="junction_aa",
    loci=("TRB",),
    use_v=True,
    use_j=True,
    use_alleles=False,
)
```

## 3. Summarize the bulk repertoire

```python
summary = iu.tl.bulk_summary(bulk)
diversity = iu.tl.alpha_diversity(bulk)
dxx = iu.tl.coverage_diversity(bulk, percentages=(20, 50, 80, 90))
clonality = iu.tl.clonality_proportion(bulk)
v_family = iu.tl.gene_usage(bulk, gene="v_call", level="family")

iu.pl.bulk_diversity(diversity, metric="shannon")
iu.pl.coverage_diversity(dxx)
iu.pl.clonality(clonality)
```

## 4. Test longitudinal expansion

```python
expansion = iu.tl.longitudinal_expansion(
    bulk,
    baseline_sample="pre",
    comparison_samples=["post"],
    fold_change=2,
)
```

The exact-test workflow is intended for initial longitudinal exploration. When
biological replicates are available, prefer the PyDESeq2 interface described in
{doc}`../user-guide/bulk-analysis`.

## 5. Link bulk clones to cells

```python
links = iu.pp.link_bulk_to_single_cell(
    bulk_chains,
    cell_chains,
    definition,
)
```

`links` retains unmatched and ambiguous mappings by default so cross-platform
disagreement remains visible instead of being silently discarded.

## 6. Continue in the single-cell ecosystem

```python
airr = iu.pp.to_scirpy(cell_chains)
mdata = iu.pp.merge_with_transcriptome(adata, airr)

iu.tl.scirpy_repertoire(mdata, groupby="sample_id")
iu.pl.clonal_expansion(mdata, groupby="sample_id")
```

See {doc}`../user-guide/single-cell` for chain QC and clonotype options, and
{doc}`../user-guide/integration` for phenotype and longitudinal analysis.
