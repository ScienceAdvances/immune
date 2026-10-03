# immune: purpose, architecture, and API

Status: implemented RNA/VDJ ecosystem foundation, October 2026.

## Purpose

immune analyzes bulk and single-cell adaptive immune receptor repertoires.
Its primary study combines single-cell RNA and VDJ data, identifies clones,
describes clone phenotypes, tracks their detection/expansion over time,
and connects MiXCR bulk clones to single-cell states.

The package incorporates PhenoTrack-style clone/state distributions and
independent-approximation flows, and CloneTrack-style longitudinal expansion.
It retains existing Python-native bulk statistics and important immunarch
features. It has no VDJtools command-line interoperability namespace.

## Ecosystem boundaries

cellscope owns RNA QC, normalization, embeddings, clustering, annotation,
functional scoring, pseudobulk, and expression inference. immune reads
their results from native objects and constructs receptor-defined groups.
The application calls cellscope directly for expression analysis. Neither
package depends on the other; native objects share annotations and graphs.

Scirpy owns mature AIRR reading, chain indexing/QC, exact clonotypes,
receptor similarity clustering, and receptor visualization. immune calls
its public methods rather than reimplementing the underlying algorithms.
Generic RNA adapters and the expression R bridge have been removed from immune.
VDJ interfaces consume existing RNA annotations for clone interpretation.

## Data model

```text
MuData (cell observations)
  mod['gex']: RNA AnnData, owned by cellscope/Scanpy
    layers['counts']: original integer counts
    obs: RNA QC, cell_type, cell_state, scores
    obsm: expression embeddings
  mod['airr']: AIRR AnnData, compatible with Scirpy
    obsm['airr']: variable-length complete receptor-chain records
    obsm['chain_indices']: Scirpy chain-selection indices
    obs: chain QC, clone identifiers, clone size/status
  obs: public sample/donor/library/condition/time metadata,
       prefixed modality annotations, has_gex and has_airr
  uns['immune']: joint result tables, parameters, backend provenance

ImmuneProject (study management / legacy bulk workflows)
  bulk_chains: canonical chain table
  single_cell_chains: canonical legacy cell-chain table
  sample_metadata: sample-level study information
  links: explicit repertoire associations
```

The bulk and single-cell axes remain distinct. Bulk rows represent measured
rearrangements; read/template/UMI counts are not fictitious cells. Legacy
canonical tables remain valid inputs. `io.read_10x` returns the existing
table format; `io.read_10x_vdj` returns native AIRR AnnData.

### Cell identity and modality alignment

Use `library_id:barcode` in both gex and airr. Keep original `barcode`,
biological `sample_id`, and `donor_id`. RNA and VDJ library IDs must refer
to the same paired cell library. Different biological samples can reuse
raw 10x barcodes but must not share composite identifiers.

Outer joins retain RNA-only and VDJ-only cells. Joint phenotype analyses use
cells with the required annotations and report their denominator. Inner
joins require `join='inner'`. Conflicting metadata for a shared cell ID
raises an error. Updating observations explicitly publishes prefixed columns
such as `gex:cell_state` and `airr:clone_id`, including with MuData >=0.4.

### Receptor and clone definitions

Full chain records are preserved in Scirpy's AIRR/Awkward representation.
Indexing chooses productive analysis chains but does not destroy the raw
records. Exact clonotypes use nucleotide sequence identity. Similarity
clusters use Scirpy's separate `define_clonotype_clusters` API.

`tl.define_clonotypes` defaults to donor-scoped, paired-arm, all-dual-chain
matching. A nonmissing donor annotation is required under this default.
Use `scope=None` deliberately for sequence sharing across donors. Scirpy
numeric cluster IDs are local to one analysis object; persist the object and
its definition together rather than interpreting IDs across independent runs.

Record sequence, metric, scope, receptor-arm/dual-chain rules and gene
constraints under `uns['immune'][key]`. Exact clone identities and bulk
single-locus match keys are different fields. A TRB can match multiple full
paired clonotypes; every association and its ambiguity remain available.

## Public namespaces and interfaces

```python
import immune as iu

iu.io       # Canonical/native repertoire readers and persistence
iu.pp       # Schema, filtering, chain QC, object assembly
iu.tl       # Repertoire, phenotype, tracking, bulk matching, joint expression
iu.pl       # Native repertoire and joint plotting
iu.get      # Result/chain/annotation extraction
iu.datasets # Deterministic paired-donor examples
```

### IO and preprocessing

| Interface | Behavior |
| --- | --- |
| `io.read_mixcr`, `io.read_airr`, `io.read_10x` | Existing canonical bulk/cell-chain readers |
| `io.read_10x_vdj(path, library_id=..., sample_id=..., donor_id=...)` | Native Scirpy CSV/JSON reader; retain raw chains |
| `io.read_airr_anndata(path, library_id=..., ...)` | Native AIRR reader with composite IDs |
| `io.read_h5ad`, `io.read_h5mu`, `io.write` | Native object persistence |
| `pp.to_scirpy(chains)` | Legacy canonical chain table to AIRR AnnData |
| `pp.merge_with_transcriptome(gex, airr, join='outer')` | Validate and assemble MuData |
| `pp.scirpy_qc(data, ...)` | Scirpy indexing and chain QC |
| `pp.CloneDefinition(...)` | Immutable single-locus key rules for bulk matching |
| `pp.clone_abundance`, `pp.filter_repertoire`, `pp.downsample_repertoire` | Legacy table aggregation/filtering/subsampling |

### Single-cell and joint tools

| Interface | Result |
| --- | --- |
| `tl.define_clonotypes(data, sequence='nt', metric='identity', scope='donor_id', ...)` | Native exact paired clonotypes or explicitly named similarity clusters |
| `tl.clone_summary(data, sample_col='sample_id', clone_key='clone_id')` | Sample/clone cell counts and frequencies |
| `tl.clonal_expansion(data, min_cells=2)` | Sample-specific clone_size and expansion-status annotations |
| `tl.phenotype_composition(data, phenotype_col='cell_state')` | Clone/state counts, within-clone fractions, sample fractions |
| `tl.phenotype_diversity(data)` | State richness, Shannon entropy, effective states |
| `tl.clone_state_enrichment(data)` | Exploratory within-sample Fisher enrichment and BH |
| `tl.phenotypic_flux(data, from_sample=..., to_sample=...)` | L1 changes in shared-clone state distributions |
| `tl.phenotype_flow(data, from_sample=..., to_sample=...)` | Independent-approximation clone-state outflow/inflow |
| `tl.track_clones(data, sample_metadata=...)` | Donor-separated count/frequency/detection trajectories |
| `tl.longitudinal_expansion(data, sample_metadata=...)` | Donor-specific baseline/follow-up tests, time-selection correction and global BH |
| `tl.match_bulk(data, bulk_chains, sample_pairs=..., definition=...)` | Cell/full-clone/single-locus/bulk association table, units and ambiguity |
| `tl.annotate_bulk_matches(data)` | Per-cell bulk annotations for a uniquely selected comparison |

Native phenotype outputs are stored in `data.uns['immune'][key]['table']`.
`params` records the denominator and analysis settings. Existing phenotype
table inputs and longitudinal table tests remain supported. Native helper
defaults use `gex` and `airr`, with explicit modality parameters where needed.

RNA expression analysis uses cellscope directly on receptor-labeled RNA
observations. immune contains no expression aggregation or inference wrappers.

### Bulk repertoire analysis

Existing APIs cover alpha/beta diversity, Hill profiles, coverage/DXX,
clonality proportions/ranks, rank-abundance, segment and V/J usage,
spectratypes, public repertoires/overlap, rarefaction, and PyDESeq2-based
sample-level differential abundance. scikit-bio supplies mature diversity
methods. Use an explicit `CloneDefinition` and count unit when comparing
pipelines or experiments.

### Plotting, extraction, datasets

Scirpy adapters remain under `pl` for repertoire networks, expansion, usage,
overlap, and receptor plots. Bulk Matplotlib plots are retained.

New `pl.phenotype_composition` produces clone/state bars;
`pl.phenotype_flow` shows a state-flow heatmap;
`pl.clone_embedding` overlays receptor annotations on the RNA embedding
without permanently changing RNA annotations. Native plotting objects are
returned; saving is explicit.

`get.chain_table` flattens full AIRR records while retaining cell metadata.
`get.airr` delegates selected indexed-chain access to Scirpy.
`get.obs_df` aligns RNA and receptor annotations.
`get.result` returns copied stored tables.

`datasets.toy_multimodal()` contains two donors, pre/post samples, paired
receptors, different TRA partners sharing TRB, and RNA-only cells.
`datasets.toy_bulk()` provides corresponding bulk rearrangements.
See `examples/rna_vdj_workflow.py` for an executable complete example.

## Statistical interpretation

- Clone sizes in single-cell summaries count cells, not receptor chains or UMIs.
- Composition fractions use clone-and-state-annotated cells as denominator.
- The phenotype-flow model estimates changes between clone-level state
  distributions; it is not observed lineage tracing of physical cells.
- Fisher state enrichment is within-sample exploratory inference. Biological
  group effects require donor/sample-level models.
- A zero/absent clone is not detected at the threshold; it does not establish
  extinction. Zero-depth samples cannot be used for expansion tests.
- Longitudinal tests run separately within donors. Selection across follow-up
  timepoints receives Bonferroni correction, then BH spans tested donor/clones.
  The inherited modified-Fisher statistic remains a CloneTrack-style method,
  not a general repeated-measures group-effect model.
- Native longitudinal sample metadata must have one sample per donor/timepoint;
  technical replicate libraries should be aggregated first.
- Bulk matching defaults to equal sample IDs. Cross-time/tissue links require
  explicit sample pairs. Available donor metadata is checked for conflicts.
- Single-locus matches can be ambiguous. Multiple bulk comparisons or chains
  per cell cannot be silently reduced to one annotation.
- Expression inference delegates raw-count aggregation and sample-level
  designs to cellscope/PyDESeq2. Biological cells are not donor replicates.

## Dependencies, sources, and extension policy

Base immune remains NumPy/Pandas/SciPy. Native receptor objects require
`immune[singlecell]`; RNA inference is installed separately through cellscope.
Plotting and other advanced backends remain optional. No commands are
launched to VDJtools or an R repertoire package.

Public backend APIs and data structures were reviewed against:

- [Scirpy receptor structure](https://scirpy.scverse.org/en/latest/data-structure.html).
- [Scirpy clone source](https://github.com/scverse/scirpy/blob/main/src/scirpy/tl/_clonotypes.py).
- [Scirpy](https://github.com/scverse/scirpy), [Scanpy](https://github.com/scverse/scanpy),
  and [MuData](https://github.com/scverse/mudata).
- [decoupler](https://github.com/saezlab/decoupler-py) and
  [PyDESeq2](https://github.com/owkin/PyDESeq2).

No third-party source files are vendored in this implementation. Future
source copying requires license review and notice retention. The local
PhenoTrack/CloneTrack reference code remains separate and unchanged.

## Future extensions

ATAC/protein modalities can be additional AnnData objects. Spatial spots
require spot-to-cell/clone mappings unless observations really are cells.
They should not be aligned by reusing unrelated barcode strings.

Antigen specificity querying, richer BCR mutation/lineage analysis, Sankey
plots, hierarchical donor models, neighborhood abundance, and large-data
optimization are future extensions rather than validated capabilities of
the initial RNA/VDJ foundation. Keep their results accessible through the
same namespaces and documented data contract.
