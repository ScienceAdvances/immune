"""Integration tests against real Scirpy, MuData, and optional cellscope."""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("anndata")
pytest.importorskip("scirpy")

import immune as iu


@pytest.fixture
def study():
    data = iu.datasets.toy_multimodal()
    iu.tl.define_clonotypes(data)
    return data


def test_native_merge_retains_unpaired_cells_and_syncs_metadata(study):
    assert study.n_obs == 96
    assert study.obs["has_airr"].sum() == 92
    assert "gex:cell_state" in study.obs
    assert "airr:clone_id" in study.obs
    assert "timepoint" in study.obs
    inner = iu.pp.merge_with_transcriptome(study.mod["gex"], study.mod["airr"], join="inner")
    assert inner.n_obs == 92


def test_conflicting_metadata_rejected(study):
    airr = study.mod["airr"].copy()
    airr.obs["donor_id"] = "other"
    with pytest.raises(ValueError, match="Conflicting donor_id"):
        iu.pp.merge_with_transcriptome(study.mod["gex"], airr)


def test_exact_pair_clonotypes_scoped_by_donor(study):
    airr = study.mod["airr"]
    d0 = set(airr.obs.loc[airr.obs["donor_id"].eq("d0"), "clone_id"])
    d1 = set(airr.obs.loc[airr.obs["donor_id"].eq("d1"), "clone_id"])
    assert not d0 & d1
    assert airr.obs.loc["lib0:cell0", "clone_id"] != airr.obs.loc["lib0:cell1", "clone_id"]
    assert airr.obs.loc["lib0:cell0", "clone_id"] == airr.obs.loc["lib1:cell0", "clone_id"]
    assert study.uns["immune"]["clone_id"]["params"]["scope"] == "donor_id"


def test_donor_scope_requires_nonmissing_labels():
    data = iu.datasets.toy_multimodal()
    del data.mod["airr"].obs["donor_id"]
    del data.obs["donor_id"]
    with pytest.raises(ValueError, match="donor_id"):
        iu.tl.define_clonotypes(data)


def test_similarity_cluster_api_uses_correct_sequence(study):
    iu.tl.define_clonotypes(study, sequence="aa", metric="identity", key_added="aa_cluster")
    assert "aa_cluster" in study.mod["airr"].obs
    assert study.uns["immune"]["aa_cluster"]["params"]["sequence"] == "aa"


def test_clone_sizes_count_cells_not_chains(study):
    table = iu.tl.clone_summary(study)
    assert table["count"].sum() == 92
    assert np.allclose(table.groupby("sample_id")["frequency"].sum(), 1)
    iu.tl.clonal_expansion(study)
    assert study.mod["airr"].obs["clone_size"].notna().all()
    assert "airr:clonal_expansion" in study.obs


def test_phenotype_composition_flux_diversity_flow(study):
    table = iu.tl.phenotype_composition(study)
    assert table["cell_count"].sum() == 92
    assert np.allclose(table.groupby(["sample_id", "clone_id"])["within_clone_fraction"].sum(), 1)
    assert iu.get.result(study, "phenotype_composition").equals(table)
    flux = iu.tl.phenotypic_flux(study, from_sample="s0", to_sample="s1")
    assert flux["phenotypic_flux"].between(0, 2).all()
    diversity = iu.tl.phenotype_diversity(study)
    assert diversity["effective_states"].between(1, 2).all()
    flow = iu.tl.phenotype_flow(study, from_sample="s0", to_sample="s1")
    assert np.isclose(flow["outflow"].sum(), 1)
    tested = iu.tl.clone_state_enrichment(study)
    assert tested["fdr"].between(0, 1).all()


def test_tracking_never_crosses_donors(study):
    tracks = iu.tl.track_clones(study)
    assert tracks.groupby("clone_id")["donor_id"].nunique().eq(1).all()
    tested = iu.tl.longitudinal_expansion(study)
    assert set(tested["donor_id"]) == {"d0", "d1"}
    assert tested["q_value"].between(0, 1).all()
    assert tested["baseline_sample"].isin(["s0", "s2"]).all()


def test_tracking_complete_zero_detection_and_order():
    abundance = pd.DataFrame(
        {"sample_id": ["pre", "post"], "clone_id": ["a", "b"], "count": [10, 10]}
    )
    metadata = pd.DataFrame(
        {"sample_id": ["post", "pre"], "donor_id": ["d", "d"], "timepoint": [2, 1]}
    )
    result = iu.tl.track_clones(abundance, metadata)
    assert len(result) == 4
    assert result["detected"].sum() == 2
    assert not result["persistent"].any()
    assert result.loc[result["clone_id"].eq("a"), "timepoint"].tolist() == [1, 2]
    metadata["timepoint"] = ["post", "pre"]
    with pytest.raises(ValueError, match="Time must"):
        iu.tl.track_clones(abundance, metadata)


def test_bulk_preserves_pair_ambiguity_and_units(study):
    links = iu.tl.match_bulk(study, iu.datasets.toy_bulk())
    assert len(links) == 92
    assert links["matched"].all()
    assert links["ambiguous"].any()
    assert set(links["bulk_count_unit"]) == {"read"}
    assert links.loc[links["ambiguous"], "n_paired_clones"].eq(2).all()
    iu.tl.annotate_bulk_matches(study)
    assert "airr:bulk_matched" in study.obs


def test_bulk_donor_mismatch_rejected(study):
    pairs = pd.DataFrame({"bulk_sample_id": ["s0"], "sc_sample_id": ["s2"]})
    with pytest.raises(ValueError, match="cross donors"):
        iu.tl.match_bulk(study, iu.datasets.toy_bulk(), sample_pairs=pairs)


def test_bulk_multi_comparison_cannot_silently_annotate(study):
    pairs = pd.DataFrame({"bulk_sample_id": ["s0", "s1"], "sc_sample_id": ["s1", "s1"]})
    iu.tl.match_bulk(study, iu.datasets.toy_bulk(), sample_pairs=pairs)
    with pytest.raises(ValueError, match="one bulk comparison"):
        iu.tl.annotate_bulk_matches(study)


def test_airr_flattening_preserves_pairs_and_original_fields(study):
    table = iu.get.chain_table(study)
    assert len(table) == 184
    assert set(table["locus"]) == {"TRA", "TRB"}
    assert table["donor_id"].notna().all()
    assert table["umi_count"].eq(3).all()
    assert iu.get.airr(study, "locus", chain="VDJ_1").dropna().eq("TRB").all()


def test_h5mu_roundtrip_with_all_joint_results(study, tmp_path):
    iu.tl.phenotype_composition(study)
    iu.tl.phenotype_diversity(study)
    iu.tl.phenotype_flow(study, from_sample="s0", to_sample="s1")
    iu.tl.track_clones(study)
    iu.tl.match_bulk(study, iu.datasets.toy_bulk())
    iu.io.write(study, tmp_path / "study.h5mu")
    restored = iu.io.read_h5mu(tmp_path / "study.h5mu")
    assert restored.n_obs == study.n_obs
    pd.testing.assert_frame_equal(
        iu.get.result(restored, "bulk_matches").reset_index(drop=True),
        iu.get.result(study, "bulk_matches").reset_index(drop=True),
        check_dtype=False,
        check_categorical=False,
    )
    assert len(iu.get.chain_table(restored)) == 184


def test_native_10x_reader_and_composite_ids(tmp_path):
    path = tmp_path / "all_contig_annotations.csv"
    pd.DataFrame(
        {
            "barcode": ["AA-1", "AA-1"],
            "chain": ["TRA", "TRB"],
            "contig_id": ["contig1", "contig2"],
            "productive": [True, True],
            "is_cell": [True, True],
            "high_confidence": [True, True],
            "cdr3_nt": ["TGTGCTTTT", "TGTACTTTT"],
            "cdr3": ["CAF", "CTF"],
            "v_gene": ["TRAV1-1", "TRBV1"],
            "j_gene": ["TRAJ1", "TRBJ1-1"],
            "d_gene": [None, "TRBD1"],
            "c_gene": ["TRAC", "TRBC1"],
            "reads": [20, 30],
            "umis": [3, 4],
        }
    ).to_csv(path, index=False)
    airr = iu.io.read_10x_vdj(path, library_id="lib", sample_id="sample", donor_id="donor")
    assert airr.obs_names.tolist() == ["lib:AA-1"]
    assert len(iu.get.chain_table(airr)) == 2
    assert airr.obs["barcode"].item() == "AA-1"


def test_joint_cellscope_aggregation_and_embedding(study):
    cs = pytest.importorskip("cellscope")
    pytest.importorskip("decoupler")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    iu.tl.clonal_expansion(study)
    pdata = iu.tl.clone_pseudobulk(study, min_cells=1, metadata_cols=["condition", "donor_id"])
    assert (
        np.asarray(pdata.X).sum()
        == study.mod["gex"].layers["counts"][study.obs["has_airr"].to_numpy()].sum()
    )
    study.mod["gex"].obsm["X_umap"] = np.random.default_rng(0).normal(size=(96, 2))
    before = study.mod["gex"].obs.columns.copy()
    assert iu.pl.clone_embedding(study) is not None
    assert study.mod["gex"].obs.columns.equals(before)
    iu.tl.phenotype_composition(study)
    iu.tl.phenotype_flow(study, from_sample="s0", to_sample="s1")
    assert iu.pl.phenotype_flow(study) is not None
    assert iu.pl.phenotype_composition(study, sample_id="s0") is not None
    assert cs.get.obs_df(study).shape[0] == 96
    plt.close("all")
