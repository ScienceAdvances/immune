"""Download-free receptor specificity example with an explicit reference."""

import pandas as pd

import immune as iu
from immune.schema import canonicalize


def main():
    chains = pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s2"],
            "donor_id": ["d1", "d1", "d2"],
            "cell_id": ["cell1", "cell2", "cell3"],
            "locus": ["TRB"] * 3,
            "productive": [True] * 3,
            "junction": ["TGTGCCAGC", "TGTGCCAGC", "TGTGCCGGT"],
            "junction_aa": ["CAS", "CAS", "CAG"],
            "v_call": ["TRBV1"] * 3,
            "j_call": ["TRBJ1-1"] * 3,
            "umi_count": [10] * 3,
        }
    )
    chains = canonicalize(chains)
    query = iu.pp.to_scirpy(chains)
    reference_chains = chains.iloc[[0]].copy()
    reference_chains["cell_id"] = "reference_cell"
    reference = iu.pp.to_scirpy(reference_chains)
    reference.obs["epitope"] = "synthetic_epitope"
    matches = iu.tl.receptor_query(
        query, reference, receptor_arms="VDJ", include_ref_cols=["epitope"]
    )
    print(matches)
    print(query.obs[["specificity_epitope"]])
    return query


if __name__ == "__main__":
    main()
