# Plotting

The `iu.pl` namespace has two behaviors:

- bulk plotting functions accept Pandas result tables produced by `iu.tl`;
- single-cell plotting functions delegate to Scirpy and accept its AnnData or
  MuData objects.

## Bulk diversity and overlap

```python
alpha = iu.tl.alpha_diversity(bulk)
distance = iu.tl.beta_diversity(bulk, metric="braycurtis")

ax = iu.pl.bulk_diversity(alpha, metric="shannon")
ax = iu.pl.repertoire_overlap(distance, cmap="mako")
```

## Clonality and rank-abundance

```python
clonality = iu.tl.clonality_proportion(bulk)
ranked = iu.tl.rank_abundance(bulk)

iu.pl.clonality(clonality)
iu.pl.rank_abundance(ranked, log_x=True, log_y=True)
```

## DXX, Hill profiles and public clonotypes

```python
dxx = iu.tl.coverage_diversity(bulk, percentages=(20, 50, 80, 90))
hill = iu.tl.hill_diversity(bulk, orders=range(6))
public = iu.tl.public_repertoire(bulk, min_samples=2)

iu.pl.coverage_diversity(dxx)
iu.pl.hill_diversity(hill)
iu.pl.public_repertoire(public, top_n=30)
```

## Gene usage and spectratypes

```python
usage = iu.tl.gene_usage(bulk, gene="v_call", level="family")
spectra = iu.tl.spectratype(bulk)

iu.pl.segment_usage(usage)
iu.pl.spectratype(spectra)
```

## Clone trajectories

```python
iu.pl.clone_trajectories(
    bulk,
    clones=("clone_a", "clone_b"),
    sample_order=("pre", "week_6", "year_1"),
)
```

## Single-cell plots

```python
iu.pl.clonal_expansion(mdata, groupby="sample_id")
iu.pl.clonotype_network(mdata, color="cell_type")
iu.pl.vdj_usage(mdata, full_combination=False)
```

These calls return the native plotting result from the installed Scirpy
version. Consult that backend's documentation for its complete keyword set.

## Publication workflow

All Python-native bulk functions return Matplotlib axes, so figures can be
customized and saved normally:

```python
ax = iu.pl.hill_diversity(hill)
ax.set_title("TCRB diversity profile")
ax.figure.tight_layout()
ax.figure.savefig("hill-diversity.svg", bbox_inches="tight")
```

Keep plotting transformations separate from statistical result tables so the
underlying numbers remain exportable and reproducible.
