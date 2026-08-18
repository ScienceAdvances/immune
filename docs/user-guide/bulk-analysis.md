# Bulk repertoire analysis

Bulk tools consume the abundance table returned by
{py:func}`immune.pp.clone_abundance`. They return Pandas DataFrames or native
objects from mature statistical backends, so results remain inspectable and
composable.

## Basic summaries and alpha diversity

```python
summary = iu.tl.bulk_summary(bulk)
alpha = iu.tl.alpha_diversity(bulk)
```

The default scikit-bio metrics are observed features, Shannon, Simpson,
inverse Simpson, Pielou evenness, Gini index and Chao1. Select an explicit set
when a study needs a smaller, preregistered panel:

```python
alpha = iu.tl.alpha_diversity(
    bulk,
    metrics=("observed_features", "shannon", "gini_index"),
)
```

Fractional MiXCR counts are rounded by default before count-based scikit-bio
estimators are called. Set `round_counts=False` only for metrics that accept
non-integer weights and document the choice.

## Hill diversity profiles

Hill numbers provide a common scale for richness, Shannon and Simpson-like
diversity across order $q$.

```python
hill = iu.tl.hill_diversity(bulk, orders=(0, 1, 2, 3, 4, 5))
iu.pl.hill_diversity(hill)
```

Low orders emphasize rare clonotypes; higher orders emphasize dominant clones.

## DXX coverage diversity

DXX is the minimum number of most abundant clonotypes needed to occupy a
requested percentage of the repertoire.

```python
dxx = iu.tl.coverage_diversity(
    bulk,
    percentages=(20, 50, 80, 90),
)
```

The output includes both the clone count (`dxx`) and its fraction of observed
richness (`dxx_fraction`).

## Clonality and rank-abundance

```python
by_proportion = iu.tl.clonality_proportion(bulk)
by_rank = iu.tl.clonality_rank(bulk)
ranked = iu.tl.rank_abundance(bulk)

iu.pl.clonality(by_proportion)
iu.pl.rank_abundance(ranked)
```

The default proportion bins are Hyperexpanded, Large, Medium, Small, Rare and
Ultra-rare. Use the annotation functions when clone-level labels are needed:

```python
annotated = iu.tl.annotate_clonality_proportion(bulk)
rank_annotated = iu.tl.annotate_clonality_rank(bulk)
```

## Gene usage and spectratypes

```python
v_family = iu.tl.gene_usage(
    bulk,
    gene="v_call",
    level="family",
    ambiguous="first",
    weighted=True,
    normalize=True,
)

v_allele = iu.tl.gene_usage(bulk, gene="v_call", level="allele")
spectra = iu.tl.spectratype(bulk, sequence_col="junction_aa")

iu.pl.segment_usage(v_family)
iu.pl.spectratype(spectra)
```

Gene levels are `family`, `segment` and `allele`. Ambiguous calls can use the
first call, be excluded, or be kept as a combined label.

## Between-sample overlap

Use beta diversity for quantitative distances:

```python
bray_curtis = iu.tl.beta_diversity(bulk, metric="braycurtis")
iu.pl.repertoire_overlap(bray_curtis)
```

Use public overlap for presence/absence similarity:

```python
jaccard = iu.tl.public_overlap(bulk, metric="jaccard")
intersection = iu.tl.public_overlap(bulk, metric="intersection")
public = iu.tl.public_repertoire(bulk, min_samples=2)
```

Jaccard and overlap-coefficient outputs are similarities, whereas scikit-bio
beta-diversity outputs are distances. Do not mix them in the same ordination
without the appropriate transformation.

## Filtering, downsampling and rarefaction

```python
filtered = iu.pp.filter_repertoire(
    bulk,
    samples=("pre", "post"),
    loci=("TRB",),
    min_count=2,
    min_frequency=1e-5,
)

sampled = iu.pp.downsample_repertoire(
    filtered,
    depth="min",
    random_state=0,
)

rarefaction = iu.tl.rarefaction(
    filtered,
    metric="observed_features",
    n_iter=50,
    seed=0,
)
```

Downsampling delegates without-replacement sampling to scikit-bio and retains
`original_count` plus `downsample_depth` for auditability.

## Longitudinal expansion without replicates

```python
expansion = iu.tl.longitudinal_expansion(
    bulk,
    baseline_sample="pre",
    comparison_samples=("week_6", "year_1"),
    fold_change=2,
    correction="bh",
)
```

This modified Fisher workflow is useful for exploratory CloneTrack-style
analysis when there is one repertoire per time point. It does not replace a
replicate-aware model.

## Replicate-aware differential abundance

With biological replicates, use PyDESeq2:

```python
sample_metadata = pd.DataFrame(
    {
        "sample_id": ["p1_pre", "p1_post", "p2_pre", "p2_post"],
        "patient": ["p1", "p1", "p2", "p2"],
        "condition": ["pre", "post", "pre", "post"],
    }
)

results, dataset, statistics = iu.tl.differential_abundance(
    bulk,
    sample_metadata,
    design="~patient + condition",
    contrast=("condition", "post", "pre"),
    min_total_count=10,
    min_samples=2,
)

statistics.plot_MA()
```

The fitted PyDESeq2 dataset and statistics objects are returned alongside the
tidy result table so native diagnostics remain available.
