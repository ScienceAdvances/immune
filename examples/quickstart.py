"""Minimal immune workflow; replace paths and sample names with your data."""

import pandas as pd

import immune as iu

project = iu.ImmuneProject()
project.add_mixcr("data/pre/clones.tsv", sample_id="pre")
project.add_mixcr("data/post/clones.tsv", sample_id="post")
project.add_10x("data/cellranger/outs", sample_id="post_sc")

definition = iu.pp.CloneDefinition(sequence="junction", loci=("TRB",))

bulk_abundance = project.bulk_abundance(definition)
links = project.build_links(definition)
expansion = project.test_expansion(
    baseline_sample="pre",
    comparison_samples=["post"],
    definition=definition,
)

# adata.obs can be used directly in place of this DataFrame.
cell_metadata = pd.DataFrame(
    {
        "sample_id": ["post_sc"],
        "cell_id": ["AAACCTGAGACAGGCT-1"],
        "cell_state": ["effector"],
    }
)
composition = project.phenotype_composition(
    cell_metadata,
    phenotype_col="cell_state",
    definition=definition,
)

print(bulk_abundance.head())
print(links.head())
print(expansion.head())
print(composition.head())
