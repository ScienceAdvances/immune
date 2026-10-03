# Mature backends

`immune` owns data normalization, explicit clonotype keys and cross-modality
coordination. Statistical and single-cell algorithms are delegated whenever a
well-maintained Python implementation exists.

| Scientific task | Backend | `immune` entry point |
|---|---|---|
| Alpha/beta diversity | scikit-bio | `iu.tl.alpha_diversity`, `beta_diversity` |
| Hill numbers | scikit-bio | `iu.tl.hill_diversity` |
| Count downsampling | scikit-bio | `iu.pp.downsample_repertoire` |
| Replicate-aware clonotype DA | PyDESeq2 | `iu.tl.differential_abundance` |
| Receptor QC and clonotyping | Scirpy | `iu.pp.scirpy_qc`, `iu.tl.scirpy_define_clonotypes` |
| Bulk plots | Matplotlib/seaborn | `iu.pl.*` |

## Lazy optional dependencies

Backends are imported only when used. If an optional dependency is missing,
the error names the appropriate extra:

```bash
python -m pip install 'immune[bulk]'
python -m pip install 'immune[singlecell]'
```

```python
iu.backend_status()
```

## Native objects remain available

Adapters return native backend objects where those objects carry important
diagnostics or fitted state. For example, differential abundance returns the
tidy results, PyDESeq2 dataset and statistics objects; Scirpy results retain their native graph and AIRR structures.

This design keeps `immune` convenient without hiding backend-specific methods,
plots, diagnostics or serialization.

## Function ownership

Small transformations that define the package's data contract—canonicalization,
clone-key construction, tidy summaries and exact matching—are implemented in
`immune`. Complex statistical models, sequence-receptor workflows and latent
representations use mature libraries.

RNA algorithms are owned by cellscope and called directly in application code.
immune retains shared object/annotation handling for receptor-focused analysis.
