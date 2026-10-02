---
html_theme.sidebar_secondary.remove: true
html_theme.sidebar_primary.remove: true
---

# immune

<div class="immune-hero">
  <img src="_static/immune-logo.svg" alt="immune" class="immune-hero-logo">
  <p class="immune-kicker">Repertoire analysis in the Python single-cell ecosystem</p>
  <h1>Bulk TCR, single-cell V(D)J, and transcriptomes in one coherent API</h1>
  <p class="immune-lead">
    Read current MiXCR, Cell Ranger and AIRR outputs; quantify repertoire structure;
    connect clonotypes to cell states; and use mature Scanpy, Scirpy and scikit-bio backends.
  </p>
  <div class="immune-hero-actions">
    <a class="immune-button primary" href="getting-started/quickstart.html">Get started</a>
    <a class="immune-button secondary" href="reference/index.html">API reference</a>
    <a class="immune-button secondary" href="https://github.com/ScienceAdvances/immune">GitHub</a>
  </div>
</div>

::::{grid} 1 2 3 3
:gutter: 3
:class-container: immune-feature-grid

:::{grid-item-card} Bulk repertoires
:class-card: immune-card
:link: user-guide/bulk-analysis
:link-type: doc

MiXCR/AIRR ingestion, richness and Hill diversity, DXX, clonality bins,
gene usage, spectratypes, overlap, rarefaction and differential abundance.

```python
alpha = iu.tl.alpha_diversity(bulk)
dxx = iu.tl.coverage_diversity(bulk)
```
:::

:::{grid-item-card} Single-cell V(D)J
:class-card: immune-card
:link: user-guide/single-cell
:link-type: doc

Current 10x V(D)J readers and thin Scirpy adapters for chain QC, clonotypes,
expansion, repertoire overlap and receptor-aware visualization.

```python
airr = iu.pp.to_scirpy(chains)
iu.tl.scirpy_repertoire(airr)
```
:::

:::{grid-item-card} Integrated studies
:class-card: immune-card
:link: user-guide/integration
:link-type: doc

Match bulk and single-cell clonotypes, analyze longitudinal expansion and
phenotype flow, and combine receptor data with Scanpy/scvi-tools/Milo workflows.

```python
links = iu.pp.link_bulk_to_single_cell(
    bulk_chains, cell_chains, definition
)
```
:::

::::

## A Scanpy-style workflow

`immune` separates input/output, preprocessing, tools and plotting into four
predictable namespaces.

::::{grid} 1 2 4 4
:gutter: 2

:::{grid-item-card} `iu.io`
Read MiXCR, Cell Ranger V(D)J and AIRR files into one canonical chain table.
:::
:::{grid-item-card} `iu.pp`
Define clonotypes, build abundance tables, filter/downsample and convert to Scirpy.
:::
:::{grid-item-card} `iu.tl`
Run repertoire statistics, longitudinal tests, phenotype and single-cell analyses.
:::
:::{grid-item-card} `iu.pl`
Plot bulk summaries or delegate receptor-aware plots to Scirpy.
:::

::::

::::{grid} 1 1 2 2
:gutter: 4

:::{grid-item}
## Install

```bash
python -m pip install immune
python -m pip install 'immune[bulk,singlecell,plot]'
```

Use Python 3.10 or newer.

{bdg-primary}`Version 0.5.0` {bdg-secondary}`MIT` {bdg-success}`Python`
:::

:::{grid-item}
## Start from real outputs

```python
import immune as iu
import pandas as pd

pre = iu.io.read_mixcr("pre.clones.tsv", sample_id="pre")
post = iu.io.read_mixcr("post.clones.tsv", sample_id="post")
chains = pd.concat([pre, post], ignore_index=True)

definition = iu.pp.CloneDefinition(sequence="junction", loci=("TRB",))
bulk = iu.pp.clone_abundance(chains, definition)
iu.pl.clonality(iu.tl.clonality_proportion(bulk))
```
:::

::::

```{toctree}
:hidden:
:maxdepth: 2

getting-started/index
user-guide/index
reference/index
ecosystem/index
development/index
```
