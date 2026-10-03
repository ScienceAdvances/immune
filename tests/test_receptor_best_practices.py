"""Real Scirpy specificity matches and AIRR alignment preservation."""

import subprocess
import sys
from importlib.util import find_spec
from pathlib import Path

import pandas as pd
import pytest

import immune as iu
from immune.schema import canonicalize


def _chains():
    return canonicalize(
        pd.DataFrame(
            {
                "sample_id": ["s1", "s1", "s2"],
                "donor_id": ["d1", "d1", "d2"],
                "cell_id": ["a", "b", "c"],
                "locus": ["TRB"] * 3,
                "productive": [True] * 3,
                "junction": ["TGTGCCAGC", "TGTGCCAGC", "TGTGCCGGT"],
                "junction_aa": ["CAS", "CAS", "CAG"],
                "v_call": ["TRBV1"] * 3,
                "j_call": ["TRBJ1-1"] * 3,
                "umi_count": [10] * 3,
                "sequence_alignment": ["AAATGT", "AAATGT", "AAATGT"],
                "germline_alignment_d_mask": ["AAATGT", "AAATGA", "AAATGT"],
            }
        )
    )


def test_interfaces_import_without_cellscope():
    for name in (
        "receptor_query",
        "bcr_clonotypes",
        "mutational_load",
        "clonotype_network",
        "clonotype_modularity",
    ):
        assert callable(getattr(iu.tl, name))


def test_vdj_package_has_no_expression_wrappers_or_dependencies():
    for name in ("advanced", "expression", "r"):
        assert not hasattr(iu, name)
        assert find_spec(f"immune.{name}") is None
    for name in (
        "scvi",
        "milo",
        "pseudobulk",
        "clone_pseudobulk",
        "clone_expression",
        "scanpy_workflow",
        "sccoda",
        "communication",
        "nichenet",
        "velocity",
        "mixscape",
    ):
        assert not hasattr(iu.tl, name)
    config = pytest.importorskip("tomllib").loads(
        (Path(__file__).parents[1] / "pyproject.toml").read_text()
    )
    extras = config["project"]["optional-dependencies"]
    assert set(extras).isdisjoint({"advanced", "joint"})
    requirements = " ".join(
        config["project"]["dependencies"]
        + [dependency for group in extras.values() for dependency in group]
    )
    for dependency in (
        "cellscope",
        "scvi-tools",
        "pertpy",
        "decoupler",
        "muon",
        "squidpy",
        "snapatac2",
    ):
        assert dependency not in requirements


def test_vdj_import_does_not_load_rna_or_other_assay_backends():
    code = """
import sys
from importlib.abc import MetaPathFinder
class Guard(MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'cellscope', 'scanpy', 'scvi', 'pertpy', 'muon', 'squidpy', 'snapatac2', 'decoupler'}:
            raise RuntimeError('Out-of-scope import: ' + fullname)
sys.meta_path.insert(0, Guard())
import immune
assert callable(immune.tl.receptor_query)
assert callable(immune.pp.to_scirpy)
"""
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


def test_exact_query_returns_annotation_in_cell_order():
    pytest.importorskip("scirpy")
    query = iu.pp.to_scirpy(_chains())
    ref_chains = _chains().iloc[[0]].copy()
    ref_chains["cell_id"] = "reference"
    reference = iu.pp.to_scirpy(ref_chains)
    reference.obs["epitope"] = "antigen1"
    table = iu.tl.receptor_query(
        query, reference, receptor_arms="VDJ", include_ref_cols=["epitope"]
    )
    assert table.index.equals(query.obs_names)
    assert query.obs.loc["a", "specificity_epitope"] == "antigen1"
    assert query.obs.loc["b", "specificity_epitope"] == "antigen1"
    assert pd.isna(query.obs.loc["c", "specificity_epitope"])
    assert iu.get.result(query, "specificity").shape == table.shape


def test_mutational_load_retains_alignment_fields():
    pytest.importorskip("scirpy")
    adata = iu.pp.to_scirpy(_chains())
    assert "sequence_alignment" in adata.obsm["airr"].fields
    iu.tl.mutational_load(adata)
    assert "mutation_count" in adata.obsm["airr"].fields
    table = iu.get.chain_table(adata)
    counts = table.set_index("cell_id")["mutation_count"]
    assert counts["a"] == 0 and counts["b"] == 1


def test_mutational_load_rejects_junction_only_data():
    pytest.importorskip("scirpy")
    adata = iu.pp.to_scirpy(
        _chains().drop(columns=["sequence_alignment", "germline_alignment_d_mask"])
    )
    with pytest.raises(ValueError, match="aligned AIRR"):
        iu.tl.mutational_load(adata)


def test_bcr_clustering_respects_donor_scope():
    pytest.importorskip("scirpy")
    chains = _chains()
    chains["locus"] = "IGH"
    chains["v_call"] = "IGHV1-2"
    chains["j_call"] = "IGHJ1"
    chains["junction"] = "TGTGCCAGC"
    chains["junction_aa"] = "CAS"
    adata = iu.pp.to_scirpy(chains)
    iu.tl.bcr_clonotypes(adata, cutoff=15, clonotype_kwargs={"receptor_arms": "VDJ"})
    assert adata.obs.loc["a", "bcr_clone_id"] == adata.obs.loc["b", "bcr_clone_id"]
    assert adata.obs.loc["a", "bcr_clone_id"] != adata.obs.loc["c", "bcr_clone_id"]
