``iu.pl``: Plotting
===================

Bulk functions accept result DataFrames. Single-cell functions forward to the
installed Scirpy plotting API.

.. currentmodule:: immune.pl

Bulk plots
----------

.. autosummary::
   :toctree: generated

   bulk_diversity
   clonality
   clone_trajectories
   coverage_diversity
   hill_diversity
   public_overlap
   public_repertoire
   rank_abundance
   segment_usage
   spectratype
   repertoire_overlap

Single-cell plots
-----------------

.. autosummary::
   :toctree: generated

   alpha_diversity
   clonal_expansion
   clonotype_imbalance
   clonotype_modularity
   clonotype_network
   embedding
   group_abundance
   vdj_usage
   clone_embedding
   phenotype_composition
   phenotype_flow
