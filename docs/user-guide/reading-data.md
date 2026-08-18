# Reading repertoire data

## MiXCR

{py:func}`immune.io.read_mixcr` recognizes common `exportClones` output fields
and MiXCR AIRR exports.

```python
mixcr = iu.io.read_mixcr(
    "sample.clones.tsv",
    sample_id="sample_01",
    productive_only=True,
)
```

A small explicit export is easier to audit than a preset-dependent table:

```bash
mixcr exportClones \
  --drop-default-fields \
  -cloneId -count -fraction \
  -vHit -jHit \
  -nFeature CDR3 -aaFeature CDR3 \
  sample.clns sample.clones.tsv
```

For UMI experiments, also export `-uniqueTagCount Molecule`. AIRR is the
preferred interchange format when it includes the fields required by the
chosen clone definition.

## Cell Ranger V(D)J

{py:func}`immune.io.read_10x` accepts a file or an `outs` directory. Supported
current outputs are:

- `filtered_contig_annotations.csv`
- `all_contig_annotations.csv`
- `airr_rearrangement.tsv`
- `clonotypes.csv`

```python
chains = iu.io.read_10x(
    "cellranger/outs",
    sample_id="patient_01_week_6",
    productive_only=True,
    high_confidence_only=True,
)
```

When a directory is supplied, cell-level filtered contigs are preferred over
AIRR, all contigs and high-level clonotypes. Use a direct file path when a
different representation is intentional.

:::{note}
`clonotypes.csv` is a high-level clone table and does not contain cell barcodes.
Use filtered contigs for cell-level RNA/TCR integration.
:::

## AIRR rearrangements

```python
airr = iu.io.read_airr(
    "airr_rearrangement.tsv",
    sample_id="sample_01",
)
```

The AIRR reader provides a stable route for data produced by tools not handled
directly. It requires enough receptor information to identify a rearrangement,
such as junction sequence or gene calls.

## Combining samples

Always assign stable sample IDs before concatenation:

```python
tables = [
    iu.io.read_mixcr(path, sample_id=sample_id)
    for sample_id, path in sample_paths.items()
]
chains = pd.concat(tables, ignore_index=True)
```

After reading, check locus, productive status, missing junctions, count fields
and duplicated source IDs before proceeding.
