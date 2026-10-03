"""Regression coverage for source abundance, cell identity and workflow result keys."""

import pandas as pd
import pytest

import immune as iu
from immune.schema import canonicalize


def _chains():
    return canonicalize(
        pd.DataFrame(
            {
                "sample_id": ["s1"] * 4,
                "source_clonotype_id": ["a", "a", "b", "c"],
                "sequence_id": ["a1", "a2", "b1", "c1"],
                "locus": ["TRB"] * 4,
                "productive": [True] * 4,
                "junction": ["TGTGCCAGC"] * 3 + ["TGTGCCGGT"],
                "junction_aa": ["CAS"] * 3 + ["CAG"],
                "cell_count": [10, 10, 5, 15],
            }
        )
    )


def test_bulk_counts_deduplicate_sources_then_sum_clones():
    result = iu.clone_abundance(_chains()).set_index("junction")
    assert result.loc["TGTGCCAGC", "count"] == 15
    assert result["count"].sum() == 30
    assert result["frequency"].eq(0.5).all()


def test_bulk_record_fallback_and_anonymous_rows():
    chains = _chains()
    chains["source_clonotype_id"] = pd.NA
    chains.loc[1, "sequence_id"] = "a1"
    assert iu.clone_abundance(chains)["count"].sum() == 30
    chains["sequence_id"] = pd.NA
    assert iu.clone_abundance(chains)["count"].sum() == 40


def test_bulk_conflicting_duplicate_source_counts_are_rejected():
    chains = _chains()
    chains.loc[1, "cell_count"] = 11
    with pytest.raises(ValueError, match="Conflicting counts"):
        iu.clone_abundance(chains)


def test_bulk_source_ids_are_namespaced_by_source():
    chains = _chains().iloc[[0, 2]].copy()
    chains["source_clonotype_id"] = "same"
    chains["source"] = ["capture1", "capture2"]
    assert iu.clone_abundance(chains)["count"].item() == 15


@pytest.mark.parametrize("library_column", [False, True])
def test_scirpy_preserves_qualified_ids_and_original_barcodes(library_column):
    pytest.importorskip("scirpy")
    chains = _chains().iloc[:2].copy()
    chains["sample_id"] = ["s1", "s2"]
    chains["cell_id"] = ["s1:AAAC", "s2:AAAC"]
    chains["barcode"] = "AAAC"
    if library_column:
        chains["library_id"] = ["s1", "s2"]
    data = iu.pp.to_scirpy(chains)
    assert data.obs_names.tolist() == chains.cell_id.tolist()
    assert data.obs.barcode.tolist() == ["AAAC", "AAAC"]
    assert data.obs_names.is_unique


def test_scirpy_rejects_final_id_collisions():
    pytest.importorskip("scirpy")
    chains = _chains().iloc[:2].copy()
    chains["sample_id"] = ["s1", "s2"]
    chains["cell_id"] = "AAAC"
    with pytest.raises(ValueError, match="Final cell IDs"):
        iu.pp.to_scirpy(chains, library_id="same")
    with pytest.raises(ValueError, match="Final cell IDs"):
        iu.pp.to_scirpy(chains, make_cell_ids_unique=False)
    assert iu.pp.to_scirpy(chains).obs_names.tolist() == ["s1:AAAC", "s2:AAAC"]


@pytest.mark.parametrize("custom_key", [None, "my_clones"])
def test_repertoire_workflow_summarizes_aa_result_not_stale_nt(custom_key):
    ir = pytest.importorskip("scirpy")
    data = iu.datasets.toy_multimodal(cells_per_sample=8)
    # A stale NT result must not be read by an AA run.
    data.mod["airr"].obs["clone_id"] = "stale_clone"
    options = {} if custom_key is None else {"key_added": custom_key}
    key = custom_key or "cc_aa_identity"
    iu.tl.scirpy_repertoire(data, groupby="sample_id", sequence="aa", clonotype_kwargs=options)
    expected = data.copy()
    ir.tl.clonal_expansion(expected, target_col=key, expanded_in="sample_id")
    pd.testing.assert_series_equal(
        data.mod["airr"].obs["clonal_expansion"],
        expected.mod["airr"].obs["clonal_expansion"],
    )
    assert key in data.uns["immune"]
    assert data.uns["immune"][key]["params"]["key_added"] == key
    assert data.mod["airr"].obs["clone_id"].eq("stale_clone").all()


def test_repertoire_workflow_aa_on_fresh_object():
    pytest.importorskip("scirpy")
    data = iu.datasets.toy_multimodal(cells_per_sample=8)
    iu.tl.scirpy_repertoire(data, groupby="sample_id", sequence="aa")
    assert "cc_aa_identity" in data.mod["airr"].obs
    assert "clonal_expansion" in data.mod["airr"].obs


def test_conflicting_clonotype_output_keys_are_rejected():
    pytest.importorskip("scirpy")
    data = iu.datasets.toy_multimodal(cells_per_sample=8)
    with pytest.raises(ValueError, match="Conflicting key_added"):
        iu.tl.define_clonotypes(data, key_added="a", clonotype_kwargs={"key_added": "b"})
