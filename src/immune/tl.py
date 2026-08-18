"""Tools: repertoire statistics, phenotype analysis and mature backend models."""

from .advanced import pseudobulk, run_milo, run_scvi, scanpy_standard_workflow
from .bulk import (
    bulk_summary,
    differential_clonotype_abundance,
    gene_usage,
    rarefaction_curve,
    segment_usage,
    skbio_alpha_diversity,
    skbio_beta_diversity,
    spectratype,
)
from .phenotype import phenotype_composition, phenotype_flow, phenotypic_flux
from .repertoire import (
    annotate_clonality_proportion,
    annotate_clonality_rank,
    clonality_proportion,
    clonality_rank,
    coverage_diversity,
    hill_diversity,
    public_overlap,
    public_repertoire,
    rank_abundance,
)
from .singlecell import (
    run_scirpy_repertoire,
    scirpy_define_clonotypes,
    scirpy_summary,
)
from .stats import repertoire_metrics, test_longitudinal_expansion, timecourse

# Concise Scanpy-style names. Explicit backend names remain available when
# analysis scripts need to make provenance visible.
alpha_diversity = skbio_alpha_diversity
beta_diversity = skbio_beta_diversity
differential_abundance = differential_clonotype_abundance
dxx = coverage_diversity
longitudinal_expansion = test_longitudinal_expansion
rarefaction = rarefaction_curve
define_clonotypes = scirpy_define_clonotypes
repertoire_summary = scirpy_summary
scirpy_repertoire = run_scirpy_repertoire
scanpy_workflow = scanpy_standard_workflow
scvi = run_scvi
milo = run_milo

__all__ = [
    "alpha_diversity",
    "annotate_clonality_proportion",
    "annotate_clonality_rank",
    "beta_diversity",
    "bulk_summary",
    "clonality_proportion",
    "clonality_rank",
    "coverage_diversity",
    "define_clonotypes",
    "differential_abundance",
    "differential_clonotype_abundance",
    "dxx",
    "gene_usage",
    "hill_diversity",
    "longitudinal_expansion",
    "milo",
    "phenotype_composition",
    "phenotype_flow",
    "phenotypic_flux",
    "pseudobulk",
    "public_overlap",
    "public_repertoire",
    "rank_abundance",
    "rarefaction",
    "rarefaction_curve",
    "repertoire_metrics",
    "repertoire_summary",
    "run_milo",
    "run_scirpy_repertoire",
    "run_scvi",
    "scanpy_standard_workflow",
    "scanpy_workflow",
    "scirpy_define_clonotypes",
    "scirpy_repertoire",
    "scirpy_summary",
    "scvi",
    "segment_usage",
    "skbio_alpha_diversity",
    "skbio_beta_diversity",
    "spectratype",
    "test_longitudinal_expansion",
    "timecourse",
]
