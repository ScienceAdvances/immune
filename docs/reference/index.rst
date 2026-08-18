API reference
=============

The public API follows the Scanpy convention. Most analysis code should use
``immune.io``, ``immune.pp``, ``immune.tl`` and ``immune.pl``. Lower-level
modules are documented for users who need explicit backend or data-model
control.

.. code-block:: python

   import immune as iu

   chains = iu.io.read_mixcr("sample.clones.tsv", sample_id="sample")
   bulk = iu.pp.clone_abundance(chains)
   summary = iu.tl.bulk_summary(bulk)

.. toctree::
   :maxdepth: 2

   io
   pp
   tl
   pl
   data-model
   project
   backends

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
