# Bulk–single-cell and longitudinal integration

The integration layer connects identical clonotype keys across data sources and
then summarizes cell phenotypes or longitudinal abundance. It does not infer
cell lineage or antigen specificity.

## Link bulk and single-cell clonotypes

```python
definition = iu.pp.CloneDefinition(
    sequence="junction",
    loci=("TRB",),
)

links = iu.pp.link_bulk_to_single_cell(
    bulk_chains,
    cell_chains,
    definition,
    include_unmatched=True,
)
```

The same definition is applied to both sources. Start with a tolerant exact key
and inspect unmatched/ambiguous cases before adding V/J constraints.

## Phenotype composition

Cell metadata can be an AnnData `.obs` table or a standalone DataFrame.

```python
composition = iu.tl.phenotype_composition(
    cell_chains,
    adata.obs,
    phenotype_col="cell_state",
    definition=definition,
    sample_col="sample_id",
)
```

The result describes phenotype proportions within each clonotype and sample.

## Phenotypic flux and flow

```python
flux = iu.tl.phenotypic_flux(
    composition,
    from_sample="week_6",
    to_sample="year_1",
)

flow = iu.tl.phenotype_flow(
    composition,
    from_sample="week_6",
    to_sample="year_1",
)
```

:::{important}
These functions compare sampled phenotype distributions. They do not observe
the same physical cell changing state, and their output should not be described
as direct lineage tracing.
:::

Use patient-aware resampling or a replicate-level model before making
population-level inferential claims.

## Study container

`ImmuneProject` is a convenience container for interactive studies:

```python
project = iu.ImmuneProject(sample_metadata=sample_metadata)
project.add_mixcr("pre.clones.tsv", sample_id="pre")
project.add_mixcr("post.clones.tsv", sample_id="post")
project.add_10x("cellranger/outs", sample_id="post_sc")

bulk = project.bulk_abundance(definition)
links = project.build_links(definition)
dxx = project.coverage_diversity(definition, percentages=(50, 90))
```

For reusable pipelines, the explicit `iu.io → iu.pp → iu.tl → iu.pl` data flow
is usually easier to test and serialize.

## Recommended provenance

Record at least:

- input tool and version (MiXCR/Cell Ranger);
- reference and germline database versions;
- selected chain, sequence and gene components of the clone definition;
- count unit and any rounding/downsampling;
- sample metadata and contrasts;
- `immune` and backend package versions;
- random seeds and filtering thresholds.
