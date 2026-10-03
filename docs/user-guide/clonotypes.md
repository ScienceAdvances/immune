# Defining clonotypes

Clonotype identity is a scientific parameter, not a file-format detail.
`immune` therefore requires an explicit {py:class}`immune.pp.CloneDefinition`
whenever source-specific IDs should be replaced by a cross-sample key.

## Sequence choice

Use nucleotide CDR3 when exact rearrangement identity matters:

```python
nt_definition = iu.pp.CloneDefinition(
    sequence="junction",
    loci=("TRB",),
)
```

Use amino-acid CDR3 when synonymous nucleotide differences should collapse:

```python
aa_definition = iu.pp.CloneDefinition(
    sequence="junction_aa",
    loci=("TRB",),
)
```

## Gene-aware definitions

```python
gene_aware = iu.pp.CloneDefinition(
    sequence="junction_aa",
    loci=("TRB",),
    use_v=True,
    use_j=True,
    use_alleles=False,
)
```

Including V/J calls increases specificity but can reduce matching across
pipelines because gene assignment algorithms and references may disagree.
Ignoring alleles is often more robust when platform interoperability is the
goal.

## Locus choice

- `(\"TRB\",)` is a practical default for many TCR beta studies.
- `(\"TRA\", \"TRB\")` treats each chain sequence as a key component but does
  not by itself reconstruct paired receptors from independent bulk chains.
- BCR clustering and mutation burden use Scirpy adapters. Germline
  reconstruction and lineage inference remain upstream tasks.

## Assign and aggregate

```python
annotated = iu.pp.add_clonotype_ids(chains, gene_aware)
abundance = iu.pp.clone_abundance(chains, gene_aware)
```

The count column is selected from available read, UMI or cell support unless
`count_col` is provided explicitly:

```python
abundance = iu.pp.clone_abundance(
    chains,
    gene_aware,
    count_col="umi_count",
)
```

Inspect `count_unit` in the returned table before downstream comparisons.

For bulk data, repeated chain rows with the same `source_clonotype_id` count
once per resulting clone; distinct source clonotypes are summed. If that ID is
missing, `sequence_id` identifies a source record. Source IDs are scoped by
sample and `source`. Conflicting counts within the same source record raise an
error. Rows without either ID are treated as independent abundance records;
remove duplicated anonymous input rows before aggregation. Single-cell
abundances still count unique cells.

## Cross-platform links

```python
links = iu.pp.link_bulk_to_single_cell(
    bulk_chains,
    cell_chains,
    nt_definition,
    include_unmatched=True,
)
```

Choose the same definition for both inputs. For MiXCR-to-Cell Ranger matching,
nucleotide CDR3 plus locus without V/J is the most tolerant exact-match
starting point. Tighten the definition only when supported by the study design.
