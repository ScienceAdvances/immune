# immune

`immune` is a Python package for integrated analysis of bulk immune
receptor sequencing and paired single-cell RNA/TCR data. It starts from the
core ideas in CloneTrack and PhenoTrack, but uses current MiXCR, Cell Ranger
and AIRR-style inputs and keeps clone-definition rules explicit.

Version `0.4.0` uses a thin-orchestration design and a Scanpy-style public API:
`io` reads data, `pp` preprocesses it, `tl` runs analyses and `pl` plots
results. The recommended import is `import immune as iu`.

The full Sphinx website lives in [`docs/`](docs/index.md) and includes
installation, data-model, bulk, single-cell, integration and API documentation.

## What works now

- MiXCR `exportClones` tables with common default or explicitly exported
  CDR3, V/J, read, UMI and cell-count fields.
- MiXCR AIRR exports.
- Current Cell Ranger V(D)J `filtered_contig_annotations.csv`,
  `all_contig_annotations.csv`, `clonotypes.csv` and
  `airr_rearrangement.tsv`.
- Canonical AIRR-like chain table with original source IDs retained.
- Configurable clonotype definitions using nucleotide or amino-acid CDR3,
  locus, V/J genes and optional alleles.
- Sample-by-clonotype abundance tables with explicit count units.
- Exact bulk-to-single-cell clonotype links, including unmatched and ambiguous
  mappings.
- CloneTrack-compatible longitudinal modified Fisher testing with fold change,
  best time point, multiple-time-point adjustment and FDR.
- Pure-Python bulk summaries, scikit-bio alpha/beta diversity（包括 Pielou、
  Gini、inverse Simpson）, V/D/J usage, CDR3 spectratype, rarefaction,
  overlap matrices and seaborn plotting.
- Immunarch-style DXX and Hill diversity profiles, rank/proportion clonality
  bins, public repertoire analysis, family/segment/allele gene usage,
  filtering and reproducible downsampling.
- Replicate-aware clonotype differential abundance through PyDESeq2, including
  paired/blocked designs expressed as formulas.
- Scirpy AIRR conversion, chain QC, clonotyping, expansion, alpha diversity,
  repertoire overlap, spectratype and its complete repertoire plotting layer.
- A conventional Scanpy workflow plus thin adapters for scvi-tools, Pertpy
  Milo differential abundance and decoupler pseudobulk aggregation.
- The original lightweight richness, Shannon, Simpson and clonality functions
  remain available for compatibility and dependency-light checks.
- Clone-aware phenotype composition, phenotypic flux and PhenoTrack-style
  probabilistic inflow/outflow.
- A small AnnData adapter that writes clonotype IDs into `adata.obs` without
  requiring Scanpy as a base dependency.

## Installation for development

```bash
cd /path/to/immune
python -m pip install -e .
```

Install only the analysis families required by a project:

```bash
python -m pip install -e '.[bulk]'
python -m pip install -e '.[differential]'
python -m pip install -e '.[singlecell,plot]'
python -m pip install -e '.[advanced]'
# or everything
python -m pip install -e '.[all]'
```

The normal analysis path is Python-only. `immune.backend_status()` reports
which optional backends are visible.

The PyPI distribution and import namespace are both named `immune`:

```python
import immune as iu
```

## MiXCR export

`immune` recognizes many MiXCR default columns. A small, predictable export
is preferable:

```bash
mixcr exportClones \
  --drop-default-fields \
  -cloneId -count -fraction \
  -vHit -jHit \
  -nFeature CDR3 -aaFeature CDR3 \
  sample.clns sample.clones.tsv
```

For UMI data, also export `-uniqueTagCount Molecule`. AIRR export is supported
directly and is the preferred interchange format when it contains the fields
needed for the intended clone definition.

## Quick start

```python
import immune as iu
import pandas as pd

pre = iu.io.read_mixcr("pre.clones.tsv", sample_id="pre")
post = iu.io.read_mixcr("post.clones.tsv", sample_id="post")
single_cell = iu.io.read_10x("cellranger/outs", sample_id="post_sc")
bulk_chains = pd.concat([pre, post], ignore_index=True)

# Cross-platform default: beta-chain nucleotide CDR3, ignoring V/J calls.
clone_def = iu.pp.CloneDefinition(sequence="junction", loci=("TRB",))

bulk = iu.pp.clone_abundance(bulk_chains, clone_def)
links = iu.pp.link_bulk_to_single_cell(bulk_chains, single_cell, clone_def)
expanded = iu.tl.longitudinal_expansion(
    bulk,
    baseline_sample="pre",
    comparison_samples=["post"],
    fold_change=2,
)

# Mature bulk statistics.
alpha = iu.tl.alpha_diversity(bulk)
bray_curtis = iu.tl.beta_diversity(bulk, metric="braycurtis")
iu.pl.bulk_diversity(alpha, metric="shannon")
```

`ImmuneProject` remains available as a convenience container, but new analysis
code should prefer the namespaced `io → pp → tl → pl` data flow.

## Python-native bulk analysis

These functions expose familiar immunarch-style results while using
Python objects throughout:

```python
summary = iu.tl.bulk_summary(bulk)
diversity = iu.tl.alpha_diversity(bulk)
distances = iu.tl.beta_diversity(bulk, metric="braycurtis")
v_usage = iu.tl.segment_usage(bulk, segment="v_call")
spectra = iu.tl.spectratype(bulk, sequence_col="junction_aa")
rarefaction = iu.tl.rarefaction(bulk, n_iter=50)

iu.pl.bulk_diversity(diversity, metric="shannon")
iu.pl.segment_usage(v_usage)
iu.pl.spectratype(spectra)
iu.pl.repertoire_overlap(distances)
```

The core immunarch-style analyses use tidy Pandas outputs:

```python
dxx = iu.tl.coverage_diversity(bulk, percentages=(20, 50, 80, 90))
hill = iu.tl.hill_diversity(bulk, orders=range(6))
prop_bins = iu.tl.clonality_proportion(bulk)
rank_bins = iu.tl.clonality_rank(bulk)
ranked = iu.tl.rank_abundance(bulk, limit=100)

v_family = iu.tl.gene_usage(bulk, gene="v_call", level="family")
public = iu.tl.public_repertoire(bulk, min_samples=2)
jaccard = iu.tl.public_overlap(bulk, metric="jaccard")

filtered = iu.pp.filter_repertoire(bulk, min_count=2, min_frequency=1e-5)
sampled = iu.pp.downsample_repertoire(filtered, depth="min", random_state=0)

iu.pl.coverage_diversity(dxx)
iu.pl.hill_diversity(hill)
iu.pl.clonality(prop_bins)
iu.pl.rank_abundance(ranked)
iu.pl.public_repertoire(public)
iu.pl.public_overlap(jaccard)
```

With biological replicates, use PyDESeq2 instead of per-clone Fisher tests:

```python
da, dds, statistics = iu.tl.differential_abundance(
    bulk,
    sample_metadata,  # one row per sample_id
    design="~patient + condition",
    contrast=["condition", "treated", "baseline"],
    min_total_count=10,
    min_samples=2,
)

# Native PyDESeq2 diagnostics remain available.
statistics.plot_MA()
```

Use stricter within-platform definitions when appropriate:

```python
strict = iu.pp.CloneDefinition(
    sequence="junction",
    loci=("TRB",),
    use_v=True,
    use_j=True,
    use_alleles=False,
)
```

## Cell phenotype analysis

`adata.obs` can be passed directly as the cell metadata table:

```python
composition = iu.tl.phenotype_composition(
    single_cell,
    adata.obs,
    phenotype_col="cell_state",
    definition=clone_def,
    sample_col="sample_id",
)

flow = iu.tl.phenotype_flow(
    composition,
    from_sample="week_6",
    to_sample="year_1",
)
```

The flow is a clone-level independent approximation between sampled state
distributions. It must not be interpreted as observation of the same physical
cell changing state.

For a full Scirpy workflow:

```python
airr = iu.pp.to_scirpy(single_cell)
mdata = iu.pp.merge_with_transcriptome(adata, airr)

iu.tl.scirpy_repertoire(mdata, groupby="sample_id")
iu.pl.clonal_expansion(mdata, groupby="sample_id")
iu.pl.repertoire_overlap(mdata)
iu.pl.vdj_usage(mdata, full_combination=False)
```

The workflow is also decomposed into `scirpy_qc`,
`scirpy_define_clonotypes` and `scirpy_summary`, so a study can change any
scientific parameter instead of being locked into a hidden pipeline.
When a barcode appears in more than one sample, `to_scirpy()` automatically
uses `sample_id:barcode`; otherwise the original barcode is preserved so it
can align directly with the transcriptome object.

## Transcriptome and advanced analysis

```python
# Optional conventional Scanpy preprocessing and embedding.
iu.tl.scanpy_workflow(adata, layer="counts")

# Batch-aware latent representation from scvi-tools.
model = iu.tl.scvi(adata, layer="counts", batch_key="sample_id")

# Neighborhood differential abundance through Pertpy Milo.
mdata, milo = iu.tl.milo(
    adata,
    sample_col="sample_id",
    design="~ condition",
    neighbors_kwargs={"use_rep": "X_scVI"},
)
milo.plot_nhood_graph(mdata, alpha=0.1)

# Sample-by-cell-type count matrices through decoupler.
pb = iu.tl.pseudobulk(
    adata,
    sample_col="sample_id",
    groups_col="cell_type",
    layer="counts",
)
```

## Important assumptions

- Count units are not interchangeable. Bulk MiXCR data commonly uses reads or
  molecules; Cell Ranger clonotypes use cells. `immune` records the chosen
  unit in abundance tables.
- MiXCR can report fractional corrected clone counts. The longitudinal exact
  test rounds them to the nearest integer only for Fisher's test.
- scikit-bio alpha diversity rounds fractional corrected counts by default for
  count-based richness estimators such as Chao1.
- The default bulk-to-single-cell match ignores V/J because pipelines may call
  the same rearrangement differently. Set `use_v` and `use_j` when a stricter
  definition is scientifically justified.
- A cell with multiple qualifying TRB chains is retained in the raw chain
  table. Cell annotation currently selects the chain with the most UMI/read
  support and does not silently count the cell twice.
- The compatibility layer uses the modified Fisher model for initial
  usability. Replicate-aware beta-binomial/mixed models are planned.

## Tests

The test suite uses the Python standard library test runner. NumPy, Pandas and
SciPy are installed as package dependencies:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

## Documentation website

```bash
python -m pip install -e '.[docs]'
python -m sphinx -W --keep-going -b html docs docs/_build/html
```

The website uses Sphinx 9.1 and PyData Sphinx Theme 0.20, with the actual build
versions displayed in the footer.

## Planned next steps

完整的 immunarch 功能逐项对照见
[`docs/immunarch_coverage.md`](docs/immunarch_coverage.md)。

1. Repertoire-level ordination、metadata-aware PERMANOVA 和 sample feature table。
2. tcrdist3/GLIPH2 and antigen-database adapters.
3. Sample-aware statistical tests for phenotype flow and clone kinetics.
4. Reproducible HTML study reports with provenance and backend versions.
