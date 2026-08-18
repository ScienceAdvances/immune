from __future__ import annotations

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

from immune import (
    abundance_matrix,
    bulk_summary,
    differential_clonotype_abundance,
    rarefaction_curve,
    scirpy_summary,
    segment_usage,
    skbio_alpha_diversity,
    skbio_beta_diversity,
    spectratype,
    to_scirpy,
)
from immune.schema import canonicalize


class FakeDistanceMatrix:
    def __init__(self, ids: list[str]) -> None:
        self.ids = ids

    def to_data_frame(self) -> pd.DataFrame:
        return pd.DataFrame(np.eye(len(self.ids)), index=self.ids, columns=self.ids)


class FakeAirrCell(dict):
    def __init__(self, cell_id: str) -> None:
        super().__init__()
        self.cell_id = cell_id
        self.chains: list[dict[str, object]] = []

    def add_chain(self, chain: dict[str, object]) -> None:
        self.chains.append(chain)


class MatureBackendTests(unittest.TestCase):
    def setUp(self) -> None:
        self.abundance = pd.DataFrame(
            {
                "sample_id": ["pre", "pre", "post", "post"],
                "clone_id": ["a", "b", "a", "c"],
                "locus": ["TRB"] * 4,
                "junction": ["AAA", "BBB", "AAA", "CCC"],
                "junction_aa": ["K", "N", "K", "P"],
                "v_call": ["TRBV1", "TRBV2", "TRBV1", "TRBV3"],
                "j_call": ["TRBJ1", "TRBJ2", "TRBJ1", "TRBJ2"],
                "count": [10.0, 5.0, 4.0, 6.0],
                "count_unit": ["read"] * 4,
                "frequency": [2 / 3, 1 / 3, 0.4, 0.6],
            }
        )

    def test_abundance_matrix(self) -> None:
        matrix = abundance_matrix(self.abundance)
        self.assertEqual(matrix.shape, (2, 3))
        self.assertEqual(float(matrix.loc["pre", "c"]), 0.0)

    def test_python_native_segment_usage_and_spectratype(self) -> None:
        usage = segment_usage(self.abundance)
        spectra = spectratype(self.abundance)
        self.assertAlmostEqual(float(usage.loc["pre"].sum()), 1.0)
        self.assertEqual(set(spectra["length"]), {1})
        totals = spectra.groupby("sample_id")["frequency"].sum()
        self.assertTrue(np.allclose(totals, 1.0))

    def test_skbio_adapters_delegate(self) -> None:
        calls: list[tuple[str, str]] = []

        def alpha(metric, counts, ids, **kwargs):
            calls.append(("alpha", metric))
            return pd.Series(np.arange(len(ids), dtype=float), index=ids)

        def beta(metric, counts, ids, **kwargs):
            calls.append(("beta", metric))
            return FakeDistanceMatrix(ids)

        fake = SimpleNamespace(alpha_diversity=alpha, beta_diversity=beta)
        with patch("immune.bulk.require_dependency", return_value=fake):
            alpha_result = skbio_alpha_diversity(
                self.abundance, metrics=("shannon", "chao1")
            )
            beta_result = skbio_beta_diversity(self.abundance)
        self.assertEqual(set(alpha_result["metric"]), {"shannon", "chao1"})
        self.assertEqual(beta_result.shape, (2, 2))
        self.assertEqual(calls, [("alpha", "shannon"), ("alpha", "chao1"), ("beta", "braycurtis")])

    def test_bulk_summary_uses_skbio_results(self) -> None:
        diversity = pd.DataFrame(
            {
                "sample_id": ["pre", "post"],
                "metric": ["shannon", "shannon"],
                "value": [0.5, 0.7],
            }
        )
        with patch("immune.bulk.skbio_alpha_diversity", return_value=diversity):
            result = bulk_summary(self.abundance, metrics=("shannon",))
        self.assertIn("shannon", result)
        self.assertEqual(set(result["n_clonotypes"]), {2})

    def test_pydeseq2_adapter_returns_native_objects(self) -> None:
        metadata = pd.DataFrame(
            {"sample_id": ["pre", "post"], "condition": ["before", "after"]}
        )

        class FakeDDS:
            def __init__(self, **kwargs):
                self.kwargs = kwargs
                self.fitted = False

            def deseq2(self):
                self.fitted = True

        class FakeStats:
            def __init__(self, dds, **kwargs):
                self.dds = dds
                self.kwargs = kwargs
                self.results_df = pd.DataFrame(
                    {"log2FoldChange": [1.5], "padj": [0.01]}, index=["a"]
                )

            def summary(self):
                return None

        modules = [
            SimpleNamespace(DeseqDataSet=FakeDDS),
            SimpleNamespace(DeseqStats=FakeStats),
        ]
        with patch("immune.bulk.require_dependency", side_effect=modules):
            results, dds, statistics = differential_clonotype_abundance(
                self.abundance,
                metadata,
                design="~condition",
                contrast=("condition", "after", "before"),
                min_total_count=1,
                min_samples=1,
            )
        self.assertTrue(dds.fitted)
        self.assertEqual(results.loc[0, "clone_id"], "a")
        self.assertIs(statistics.dds, dds)

    def test_rarefaction_delegates_subsampling_and_metric(self) -> None:
        def subsample_counts(counts, n, replace, seed):
            result = np.zeros_like(counts)
            result[0] = n
            return result

        def alpha_diversity(metric, counts):
            return np.asarray([(np.asarray(counts) > 0).sum()], dtype=float)

        modules = [
            SimpleNamespace(alpha_diversity=alpha_diversity),
            SimpleNamespace(subsample_counts=subsample_counts),
        ]
        with patch("immune.bulk.require_dependency", side_effect=modules):
            result = rarefaction_curve(
                self.abundance,
                depths=(5,),
                n_iter=3,
            )
        self.assertEqual(len(result), 2)
        self.assertTrue((result["mean"] == 1.0).all())

    def test_scirpy_summary_delegates_current_api(self) -> None:
        calls: list[tuple[str, dict[str, object]]] = []

        def record(name):
            def inner(*args, **kwargs):
                calls.append((name, kwargs))

            return inner

        fake = SimpleNamespace(
            tl=SimpleNamespace(
                clonal_expansion=record("clonal_expansion"),
                alpha_diversity=record("alpha_diversity"),
                repertoire_overlap=record("repertoire_overlap"),
                spectratype=record("spectratype"),
            )
        )
        data = SimpleNamespace(uns={})
        with patch("immune.singlecell.require_dependency", return_value=fake):
            result = scirpy_summary(data, groupby="sample_id")
        self.assertIs(result, data)
        spectratype_call = next(kwargs for name, kwargs in calls if name == "spectratype")
        self.assertEqual(spectratype_call["target_col"], "sample_id")
        self.assertIn("spectratype_sample_id", data.uns)

    def test_to_scirpy_preserves_chains_and_disambiguates_barcodes(self) -> None:
        chains = canonicalize(
            pd.DataFrame(
                {
                    "sample_id": ["s1", "s2"],
                    "cell_id": ["AA-1", "AA-1"],
                    "sequence_id": ["c1", "c2"],
                    "source": ["10x_contig", "10x_contig"],
                    "locus": ["TRB", "TRB"],
                    "productive": [True, True],
                    "junction": ["AAA", "BBB"],
                    "junction_aa": ["K", "N"],
                    "v_call": ["TRBV1", "TRBV2"],
                    "j_call": ["TRBJ1", "TRBJ2"],
                    "read_count": [20, 30],
                    "umi_count": [2, 3],
                }
            )
        )
        captured: list[FakeAirrCell] = []

        def from_airr_cells(cells):
            captured.extend(cells)
            return "airr-adata"

        fake = SimpleNamespace(
            io=SimpleNamespace(AirrCell=FakeAirrCell, from_airr_cells=from_airr_cells)
        )
        with patch("immune.singlecell.require_dependency", return_value=fake):
            result = to_scirpy(chains)
        self.assertEqual(result, "airr-adata")
        self.assertEqual([cell.cell_id for cell in captured], ["s1:AA-1", "s2:AA-1"])
        self.assertEqual(captured[0].chains[0]["duplicate_count"], 2)


if __name__ == "__main__":
    unittest.main()
