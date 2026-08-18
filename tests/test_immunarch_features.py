from __future__ import annotations

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

import immune as iu
from immune.bulk import DEFAULT_ALPHA_METRICS


class ImmunarchFeatureTests(unittest.TestCase):
    def setUp(self) -> None:
        self.abundance = pd.DataFrame(
            {
                "sample_id": ["s1", "s1", "s1", "s2", "s2", "s3", "s3"],
                "clone_id": ["a", "b", "c", "a", "d", "a", "b"],
                "locus": ["TRB"] * 7,
                "junction": ["AAA", "BBB", "CCC", "AAA", "DDD", "AAA", "BBB"],
                "junction_aa": ["K", "N", "P", "K", "D", "K", "N"],
                "v_call": [
                    "TRBV7-9*01",
                    "TRBV7-8*02",
                    "TRBV5-1*01",
                    "TRBV7-9*01",
                    "TRBV6-5*01",
                    "TRBV7-9*01",
                    "TRBV7-8*02",
                ],
                "d_call": ["TRBD1*01"] * 7,
                "j_call": ["TRBJ2-7*01"] * 7,
                "c_call": ["TRBC2*01"] * 7,
                "count": [60, 30, 10, 50, 50, 1, 99],
                "frequency": [0.6, 0.3, 0.1, 0.5, 0.5, 0.01, 0.99],
                "count_unit": ["read"] * 7,
            }
        )

    def test_dxx_and_hill_profile(self) -> None:
        self.assertTrue(
            {"pielou_e", "gini_index", "inv_simpson"}.issubset(DEFAULT_ALPHA_METRICS)
        )
        coverage = iu.tl.coverage_diversity(self.abundance, percentages=(50, 90))
        s1 = coverage[coverage["sample_id"] == "s1"].set_index("percentage")
        self.assertEqual(int(s1.loc[50.0, "dxx"]), 1)
        self.assertEqual(int(s1.loc[90.0, "dxx"]), 2)

        def alpha(metric, counts, ids, order):
            self.assertEqual(metric, "hill")
            return pd.Series([order + 1] * len(ids), index=ids, dtype=float)

        fake = SimpleNamespace(alpha_diversity=alpha)
        with patch("immune.repertoire.require_dependency", return_value=fake):
            profile = iu.tl.hill_diversity(self.abundance, orders=(0, 1, 2))
        self.assertEqual(set(profile["q"]), {0.0, 1.0, 2.0})
        self.assertEqual(len(profile), 9)

    def test_clonality_annotations_and_summaries(self) -> None:
        bins = {"High": 0.5, "Medium": 0.1}
        annotated = iu.tl.annotate_clonality_proportion(self.abundance, bins=bins)
        self.assertEqual(
            str(annotated.query("sample_id == 's1' and clone_id == 'a'").iloc[0]["clonal_prop_bin"]),
            "High",
        )
        summary = iu.tl.clonality_proportion(self.abundance, bins=bins)
        occupied = summary.groupby("sample_id", observed=True)["occupied_frequency"].sum()
        self.assertTrue(np.allclose(occupied, 1.0))

        ranked = iu.tl.annotate_clonality_rank(self.abundance, bins=(1, 2))
        s1 = ranked[ranked["sample_id"] == "s1"].sort_values("rank")
        self.assertEqual(s1["clonal_rank_bin"].astype(str).tolist(), ["1-1", "2-2", ">2"])

    def test_public_repertoire_and_overlap(self) -> None:
        public = iu.tl.public_repertoire(self.abundance, min_samples=2)
        self.assertEqual(set(public["clone_id"]), {"a", "b"})
        self.assertEqual(int(public.set_index("clone_id").loc["a", "n_samples"]), 3)

        intersection = iu.tl.public_overlap(self.abundance, metric="intersection")
        self.assertEqual(float(intersection.loc["s1", "s2"]), 1.0)
        self.assertEqual(float(intersection.loc["s1", "s3"]), 2.0)
        jaccard = iu.tl.public_overlap(self.abundance, metric="jaccard")
        self.assertAlmostEqual(float(jaccard.loc["s1", "s2"]), 0.25)

    def test_gene_usage_levels(self) -> None:
        family = iu.tl.gene_usage(self.abundance, gene="v_call", level="family")
        self.assertIn("TRBV7", family.columns)
        self.assertAlmostEqual(float(family.loc["s1"].sum()), 1.0)
        allele = iu.tl.gene_usage(self.abundance, gene="v_call", level="allele")
        self.assertIn("TRBV7-9*01", allele.columns)

    def test_filter_and_downsample(self) -> None:
        filtered = iu.pp.filter_repertoire(
            self.abundance,
            samples=("s1",),
            min_frequency=0.2,
            min_length=1,
            sequence_col="junction_aa",
        )
        self.assertEqual(set(filtered["clone_id"]), {"a", "b"})

        def subsample(counts, n, replace, seed):
            result = np.zeros_like(counts)
            remaining = n
            for index, count in enumerate(counts):
                selected = min(int(count), remaining)
                result[index] = selected
                remaining -= selected
            return result

        fake = SimpleNamespace(subsample_counts=subsample)
        with patch("immune.preprocessing.require_dependency", return_value=fake):
            sampled = iu.pp.downsample_repertoire(self.abundance, depth=50)
        totals = sampled.groupby("sample_id")["count"].sum()
        frequencies = sampled.groupby("sample_id")["frequency"].sum()
        self.assertTrue((totals == 50).all())
        self.assertTrue(np.allclose(frequencies, 1.0))

    def test_scanpy_style_plot_api_is_exposed(self) -> None:
        for name in (
            "clonality",
            "coverage_diversity",
            "hill_diversity",
            "public_overlap",
            "public_repertoire",
            "rank_abundance",
        ):
            self.assertTrue(callable(getattr(iu.pl, name)))


if __name__ == "__main__":
    unittest.main()
