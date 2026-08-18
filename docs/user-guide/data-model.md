# Data model

`immune` uses two deliberately separate table types:

1. a **chain table** with one rearrangement/contig observation per row; and
2. an **abundance table** with one sample–clonotype observation per row.

Keeping these tables separate prevents cell-level observations, reads,
molecules and clone counts from being treated as interchangeable units.

## Canonical chain table

Every reader normalizes its output into an AIRR-like Pandas DataFrame. Important
columns are:

| Column | Meaning |
|---|---|
| `sample_id` | User-supplied or source-derived sample identifier |
| `cell_id` | Cell barcode for cell-level receptor data |
| `sequence_id` | Source contig, sequence or clone record identifier |
| `source` | Reader/source label such as `mixcr`, `10x_contig` or `airr` |
| `source_clonotype_id` | Original pipeline clonotype identifier |
| `locus` | TRA, TRB, TRG, TRD, IGH, IGK or IGL |
| `productive` | Productive rearrangement flag when available |
| `high_confidence` | Source confidence flag when available |
| `junction` | CDR3/junction nucleotide sequence |
| `junction_aa` | CDR3/junction amino-acid sequence |
| `v_call`, `d_call`, `j_call`, `c_call` | Normalized gene calls |
| `read_count` | Read support |
| `umi_count` | Molecule/UMI support |
| `cell_count` | Number of cells represented by a high-level clone row |
| `frequency` | Source-provided frequency when available |

Missing source fields remain missing; they are not invented. The raw source ID
is retained so results can be traced back to MiXCR or Cell Ranger.

```python
from immune.schema import validate_rearrangements

validate_rearrangements(chains)
```

## Abundance table

{py:func}`immune.pp.clone_abundance` aggregates the chain table using a
{py:class}`immune.pp.CloneDefinition`. The result contains at least:

| Column | Meaning |
|---|---|
| `sample_id` | Repertoire/sample identifier |
| `clone_id` | Stable hash derived from the selected clone key |
| `count` | Selected count unit aggregated for the clone |
| `frequency` | Within-sample normalized abundance |
| `count_unit` | `read`, `umi`, `cell` or another explicitly selected unit |

Sequence and gene columns used in the definition are retained for reporting.

## Count-unit rules

:::{warning}
Reads, UMIs and cells answer different questions. Do not concatenate or compare
their raw counts as if they shared the same sampling process.
:::

- Bulk MiXCR data commonly measures reads or molecules.
- Cell Ranger clonotype tables represent cells.
- Corrected MiXCR clone counts may be fractional.
- Count-based richness estimators and exact tests may require integer counts;
  functions document when rounding occurs.

Use frequencies for descriptive cross-platform comparisons, or a study design
that explicitly models each count-generating process.

## AnnData and MuData

`immune` does not replace the scverse object model. Cell-level receptor chains
can be converted to Scirpy AIRR AnnData and merged with transcriptomes into
MuData. Transcriptomic matrices, embeddings and cell annotations remain in the
native Scanpy/Scirpy objects.
