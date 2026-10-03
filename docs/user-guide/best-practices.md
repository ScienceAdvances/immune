# Single-cell Best Practices and RNA/VDJ analysis

immune 0.6 adds receptor specificity queries, BCR sequence-similarity
clustering, mutation burden and clone graph analysis for the adaptive immune
chapters of [Single-cell Best Practices](https://www.sc-best-practices.org/).
RNA analysis lives in cellscope 1.3 and works on the same AnnData/MuData
objects. immune implements V(D)J repertoire analysis and uses Python Scirpy.
It reads existing RNA annotations and graphs for receptor-focused interpretation.

## Installation

```bash
python -m pip install -e '/Users/tim/Documents/Repo/immune[best-practices,plot]'
```

`immune[best-practices]` installs the AIRR object and Scirpy backends. Importing
immune does not require cellscope or R. `iu.backend_status()` inspects the
supported V(D)J backends. Install cellscope separately for RNA analysis.

## Receptor interfaces

| Tutorial task | Interface | Backend and result |
|---|---|---|
| Read receptor data | `io.read_airr`, `io.read_10x`, `pp.to_scirpy` | Named AIRR chains, Scirpy AnnData |
| Chain QC | `pp.scirpy_qc` | Scirpy chain-pairing annotations |
| TCR clonotypes | `tl.scirpy_define_clonotypes` | Exact or distance-based native Scirpy groups |
| BCR similarity groups | `tl.bcr_clonotypes` | Normalized nucleotide Hamming distances and donor-scoped groups |
| Antigen/database query | `tl.receptor_query` | Explicit reference; ordered annotation table and prefixed obs columns |
| Somatic mutation load | `tl.mutational_load` | Scirpy AIRR mutation fields; requires germline alignments |
| Clone graph and layout | `tl.clonotype_network` | Native Scirpy network for `pl.clonotype_network` |
| Clone/transcriptome localization | `tl.clonotype_modularity` | Native Scirpy graph modularity result |
| Diversity, expansion and overlap | Existing `tl.scirpy_summary`, `tl.scirpy_repertoire`, diversity and expansion functions | Scirpy and Python repertoire tables |

`to_scirpy` now retains alignment and region-coordinate fields from AIRR input,
instead of reducing every chain to junction and V/J calls. Canonical metadata
and sample-aware cell identifiers are retained. Raw aligned AIRR fields are not
translated or regenerated. A 10x junction-only table cannot supply somatic
mutation estimates; `mutational_load` reports the missing required field.

## Query an explicit specificity reference

```python
import immune as iu

query_chains = iu.io.read_airr("query.airr.tsv", sample_id="library1")
reference_chains = iu.io.read_airr("reference.airr.tsv", sample_id="reference")
query = iu.pp.to_scirpy(query_chains)
reference = iu.pp.to_scirpy(reference_chains)
# Attach harmonized reference annotations by reference cell identifier.
reference.obs["epitope"] = reference_annotations["epitope"].reindex(reference.obs_names)

matches = iu.tl.receptor_query(
    query, reference, sequence="aa", metric="identity", receptor_arms="VDJ",
    include_ref_cols=["epitope"], strategy="unique-only",
)
print(query.obs["specificity_epitope"])
```

The reference may come from VDJdb, IEDB or a study-specific source after
harmonization. No external database is downloaded automatically. Exact CDR3
matches or similarity groups are specificity hypotheses, not experimental
proof; sequence conventions, V genes, paired chains and HLA affect the match.
`strategy="unique-only"` leaves conflicting annotations unresolved;
`strategy="json"` retains ambiguous reference values explicitly.

## BCR clustering and mutation burden

```python
bcr = iu.pp.to_scirpy(iu.io.read_airr("aligned_bcr.airr.tsv", sample_id="library1"))
# Supply donor_id in cell metadata before donor-scoped grouping.
iu.tl.bcr_clonotypes(
    bcr, scope="donor_id", cutoff=15,
    clonotype_kwargs={"same_v_gene": True, "receptor_arms": "VDJ"},
)
iu.tl.mutational_load(bcr, regions=["full"])
chains = iu.get.chain_table(bcr)
```

The normalized-Hamming cutoff uses Scirpy's native distance scale; it is not
a count of 15 nucleotide mismatches. Set V-gene, heavy/light-chain and donor
constraints according to the study. Use `receptor_arms="VDJ"` for heavy-chain
only data; the default paired-chain analysis expects suitable VJ chains.
The default clone key is `bcr_clone_id`, so a TCR clone annotation is not
overwritten. Mutation calculation requires `sequence_alignment` and
`germline_alignment_d_mask`; alternative AIRR keys can be supplied explicitly.
Region-specific estimates also require alignment-coordinate fields.

## Joint expression and immune analysis

```python
import cellscope as cs
import immune as iu

mdata = iu.pp.merge_with_transcriptome(rna, query)
cs.pp.normalize(mdata, mod="gex")
cs.tl.pca(mdata, mod="gex", n_comps=30)
cs.pp.neighbors(mdata, mod="gex")
cs.tl.umap(mdata, mod="gex")
iu.tl.scirpy_define_clonotypes(mdata, scope="donor_id")
iu.tl.clonotype_network(mdata)
iu.pl.clonotype_network(mdata)
iu.tl.clonotype_modularity(mdata, connectivity_key="gex:connectivities")

```

The application calls cellscope directly for RNA QC, normalization, integration,
expression inference, RNA velocity, communication and composition analysis.
immune exposes receptor processing and summaries. It has no generic RNA wrapper
namespace or R bridge. AnnData/MuData conversion and clone/phenotype associations
remain available to share existing results between packages.

## Verification

Real Scirpy tests check ordered reference matches, ambiguity/unmatched output,
preserved alignments, observed mutation counts and donor-restricted BCR groups.
Existing RNA/VDJ persistence and integration tests remain part of the regression
suite. Run `python examples/receptor_best_practices.py` for a small explicit
reference query without downloading data. For validation of RNA adapters,
consult the cellscope guide; availability does not imply
every model has been tested on a real cohort.
