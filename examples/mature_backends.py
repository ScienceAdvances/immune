"""Mature-backend workflow; replace paths and metadata columns for your study."""

import scanpy as sc

import immune as iu

project = iu.ImmuneProject()
project.add_mixcr("data/pre.clones.tsv", sample_id="pre")
project.add_mixcr("data/post.clones.tsv", sample_id="post")
project.add_10x("data/cellranger/outs", sample_id="post")

definition = iu.pp.CloneDefinition(sequence="junction", loci=("TRB",))

# Bulk: pure-Python summary tables and plots.
alpha = project.bulk_diversity(definition)
distances = project.bulk_distances(definition, metric="braycurtis")
summary = project.bulk_summary(definition)
v_usage = project.segment_usage(definition, segment="v_call")
spectra = project.spectratype(definition)
rarefaction = project.rarefaction(definition, n_iter=50)

iu.pl.bulk_diversity(alpha, metric="shannon")
iu.pl.segment_usage(v_usage)
iu.pl.spectratype(spectra)
iu.pl.repertoire_overlap(distances)

# Single cell: merge Scanpy gene expression with Scirpy AIRR data.
rna = sc.read_h5ad("data/rna.h5ad")
airr = project.to_scirpy()
mdata = iu.pp.merge_with_transcriptome(rna, airr)
iu.tl.scirpy_repertoire(mdata, groupby="sample_id")
iu.pl.clonal_expansion(mdata, groupby="sample_id")
iu.pl.vdj_usage(mdata, full_combination=False)

# Advanced analyses remain explicit and return their native objects.
scvi_model = iu.tl.scvi(rna, layer="counts", batch_key="sample_id")
milo_data, milo = iu.tl.milo(rna, sample_col="sample_id", design="~condition")
pseudobulk = iu.tl.pseudobulk(
    rna,
    sample_col="sample_id",
    groups_col="cell_type",
    layer="counts",
)

print(alpha)
print(distances)
