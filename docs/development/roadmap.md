# Roadmap

The next additions should extend scientific coverage without rebuilding mature
algorithms locally.

## Repertoire-level cohort statistics

- metadata-aware sample feature tables;
- PCoA/PCA and hierarchical clustering;
- PERMANOVA/ANOSIM for prespecified sample groups;
- public-clonotype combination statistics;
- additional validated overlap measures.

## Sequence similarity and annotation

- tcrdist3 adapters for bulk sequence distances;
- GLIPH2-compatible workflow integration;
- versioned VDJdb and McPAS annotation;
- k-mer, sequence logo and amino-acid property backends.

## BCR workflows

Connect Dandelion or Immcantation/Change-O for germline assignment, somatic
hypermutation and lineage analysis. These algorithms should not be recreated
inside `immune`.

## Reproducibility

- study-level provenance and backend version capture;
- exportable analysis manifests;
- reproducible HTML study reports;
- richer example datasets and end-to-end tutorials.
