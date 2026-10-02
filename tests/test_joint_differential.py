"""Native joint expression methods delegate to cellscope sample-level backends."""

import pytest

pytest.importorskip("scirpy")
pytest.importorskip("cellscope")
pytest.importorskip("decoupler")

import immune as iu


@pytest.mark.parametrize("method", ["pylimma", "edgepython_ql"])
def test_joint_expression_method_and_provenance(method):
    pytest.importorskip("pylimma" if method == "pylimma" else "edgepython")
    data = iu.datasets.toy_multimodal()
    iu.tl.define_clonotypes(data)
    iu.tl.clonal_expansion(data)
    table, pdata, models = iu.tl.clone_expression(
        data,
        method=method,
        design="~ donor_id + condition",
        contrast=("condition", "post", "pre"),
        min_cells=1,
        metadata_cols=["condition", "donor_id"],
    )
    assert not table.empty and len(models) > 0
    assert set(table["method"]) == {method}
    assert len(pdata) >= 4
    entry = data.uns["immune"]["clone_expression"]
    assert entry["params"]["method"] == method
    assert entry["backend"] == ("edgepython" if method.startswith("edgepython") else method)
