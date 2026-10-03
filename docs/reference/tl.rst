``iu.tl``: Tools
================

The tools namespace combines Python-native repertoire tables with thin adapters
to mature statistics and single-cell libraries.

.. currentmodule:: immune.tl

Bulk summaries and diversity
----------------------------

.. autosummary::
   :toctree: generated

   bulk_summary
   alpha_diversity
   beta_diversity
   coverage_diversity
   dxx
   hill_diversity
   repertoire_metrics
   rarefaction
   rarefaction_curve

Clonality, genes and public repertoires
----------------------------------------

.. autosummary::
   :toctree: generated

   annotate_clonality_proportion
   annotate_clonality_rank
   clonality_proportion
   clonality_rank
   rank_abundance
   gene_usage
   segment_usage
   spectratype
   public_repertoire
   public_overlap

Longitudinal and differential analysis
--------------------------------------

.. autosummary::
   :toctree: generated

   longitudinal_expansion
   test_longitudinal_expansion
   timecourse
   track_clones
   differential_abundance
   differential_clonotype_abundance

Phenotype analysis
------------------

.. autosummary::
   :toctree: generated

   phenotype_composition
   phenotype_flow
   phenotypic_flux
   phenotype_diversity
   clone_state_enrichment

Native RNA/VDJ joint analysis
-----------------------------

.. autosummary::
   :toctree: generated

   define_clonotypes
   clone_summary
   clonal_expansion
   match_bulk
   annotate_bulk_matches

Single-cell receptor backends
-----------------------------

.. autosummary::
   :toctree: generated

   scirpy_define_clonotypes
   scirpy_summary
   scirpy_repertoire

Best-practice receptor analysis
-------------------------------

.. autosummary::
   :toctree: generated

   receptor_query
   bcr_clonotypes
   mutational_load
   clonotype_network
   clonotype_modularity
