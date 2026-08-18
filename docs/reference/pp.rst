``iu.pp``: Preprocessing
========================

Preprocessing functions make clonotype identity and transformations explicit.

.. currentmodule:: immune.pp

Clonotypes and matrices
-----------------------

.. autosummary::
   :toctree: generated

   CloneDefinition
   add_clonotype_ids
   clone_abundance
   link_bulk_to_single_cell
   abundance_matrix

Filtering and downsampling
--------------------------

.. autosummary::
   :toctree: generated

   filter_repertoire
   downsample_repertoire

Single-cell conversion
----------------------

.. autosummary::
   :toctree: generated

   to_scirpy
   merge_with_transcriptome
   scirpy_qc
   add_clonotypes_to_anndata
