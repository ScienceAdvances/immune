"""Executable paired-donor RNA/VDJ/bulk workflow on deterministic toy data."""

import cellscope as cs

import immune as iu


def main():
    data = iu.datasets.toy_multimodal()
    cs.pp.qc(data, min_genes=5)
    cs.tl.workflow(data, n_top_genes=30, n_pcs=10, random_state=0)
    iu.tl.define_clonotypes(data, scope="donor_id")
    iu.tl.clonal_expansion(data)
    iu.tl.phenotype_composition(data)
    iu.tl.phenotype_diversity(data)
    iu.tl.clone_state_enrichment(data)
    iu.tl.phenotypic_flux(data, from_sample="s0", to_sample="s1")
    iu.tl.phenotype_flow(data, from_sample="s0", to_sample="s1")
    iu.tl.track_clones(data)
    iu.tl.longitudinal_expansion(data)
    iu.tl.match_bulk(data, iu.datasets.toy_bulk())
    iu.tl.annotate_bulk_matches(data)
    iu.pl.clone_embedding(data)
    iu.pl.phenotype_flow(data)
    print(iu.get.result(data, "clone_summary"))
    return data


if __name__ == "__main__":
    main()
