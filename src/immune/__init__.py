"""immune: Scanpy-style bulk and single-cell immune repertoire analysis."""

from . import io, pl, pp, tl
from ._optional import backend_status
from .advanced import pseudobulk, run_milo, run_scvi, scanpy_standard_workflow
from .bulk import (
    abundance_matrix,
    bulk_summary,
    differential_clonotype_abundance,
    gene_usage,
    rarefaction_curve,
    segment_usage,
    skbio_alpha_diversity,
    skbio_beta_diversity,
    spectratype,
)
from .clonotypes import (
    CloneDefinition,
    add_clonotype_ids,
    clone_abundance,
    link_bulk_to_single_cell,
)
from .io import read_10x, read_airr, read_mixcr
from .phenotype import phenotype_composition, phenotype_flow, phenotypic_flux
from .plotting import (
    plot_bulk_diversity,
    plot_clone_trajectories,
    plot_repertoire_overlap,
    plot_segment_usage,
    plot_spectratype,
)
from .preprocessing import downsample_repertoire, filter_repertoire
from .project import ImmuneProject
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
    merge_with_transcriptome,
    plot_scirpy,
    run_scirpy_repertoire,
    scirpy_define_clonotypes,
    scirpy_qc,
    scirpy_summary,
    to_scirpy,
)
from .stats import repertoire_metrics, test_longitudinal_expansion, timecourse

__all__ = [
    "CloneDefinition",
    "ImmuneProject",
    "abundance_matrix",
    "add_clonotype_ids",
    "annotate_clonality_proportion",
    "annotate_clonality_rank",
    "backend_status",
    "bulk_summary",
    "clonality_proportion",
    "clonality_rank",
    "clone_abundance",
    "coverage_diversity",
    "differential_clonotype_abundance",
    "downsample_repertoire",
    "filter_repertoire",
    "gene_usage",
    "hill_diversity",
    "io",
    "link_bulk_to_single_cell",
    "merge_with_transcriptome",
    "phenotype_composition",
    "phenotype_flow",
    "phenotypic_flux",
    "pl",
    "plot_bulk_diversity",
    "plot_clone_trajectories",
    "plot_repertoire_overlap",
    "plot_scirpy",
    "plot_segment_usage",
    "plot_spectratype",
    "pp",
    "pseudobulk",
    "public_overlap",
    "public_repertoire",
    "rank_abundance",
    "rarefaction_curve",
    "read_10x",
    "read_airr",
    "read_mixcr",
    "repertoire_metrics",
    "run_milo",
    "run_scirpy_repertoire",
    "run_scvi",
    "scanpy_standard_workflow",
    "scirpy_define_clonotypes",
    "scirpy_qc",
    "scirpy_summary",
    "segment_usage",
    "skbio_alpha_diversity",
    "skbio_beta_diversity",
    "spectratype",
    "test_longitudinal_expansion",
    "timecourse",
    "tl",
    "to_scirpy",
]

__version__ = "0.4.0"
