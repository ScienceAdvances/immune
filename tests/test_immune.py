from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from immune import (
    CloneDefinition,
    ImmuneProject,
    clone_abundance,
    link_bulk_to_single_cell,
    phenotype_composition,
    phenotype_flow,
    phenotypic_flux,
    read_10x,
    read_airr,
    read_mixcr,
    test_longitudinal_expansion,
)


class ReaderTests(unittest.TestCase):
    def test_mixcr_default_export(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "clones.tsv"
            pd.DataFrame(
                {
                    "cloneId": [1, 2, 3],
                    "cloneCount": [100.0, 25.0, 10.0],
                    "cloneFraction": [0.74, 0.19, 0.07],
                    "allVHitsWithScore": [
                        "TRBV7-9*01(400.2)",
                        "TRBV5-1*01(300.1)",
                        "TRBV6-5*01(200.0)",
                    ],
                    "allJHitsWithScore": [
                        "TRBJ2-7*01(100.0)",
                        "TRBJ1-2*01(99.0)",
                        "TRBJ2-1*01(80.0)",
                    ],
                    "nSeqCDR3": ["TGTGCCAGC", "TGTGCCAAA", "TGTGCCGGG"],
                    "aaSeqCDR3": ["CASSF", "CASKF", "CAS*F"],
                }
            ).to_csv(path, sep="\t", index=False)
            result = read_mixcr(path, sample_id="bulk_pre")
        self.assertEqual(len(result), 2)
        self.assertEqual(result.loc[0, "v_call"], "TRBV7-9*01")
        self.assertEqual(result.loc[0, "locus"], "TRB")
        self.assertEqual(result.loc[0, "sample_id"], "bulk_pre")

    def test_10x_filtered_contigs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "filtered_contig_annotations.csv"
            pd.DataFrame(
                {
                    "sample": ["sc1", "sc1", "sc1"],
                    "barcode": ["AA-1", "BB-1", "CC-1"],
                    "is_cell": [True, True, True],
                    "contig_id": ["AA_contig_1", "BB_contig_1", "CC_contig_1"],
                    "high_confidence": [True, True, False],
                    "chain": ["TRB", "TRB", "TRB"],
                    "v_gene": ["TRBV7-9", "TRBV5-1", "TRBV6-5"],
                    "d_gene": ["TRBD1", "TRBD1", "TRBD2"],
                    "j_gene": ["TRBJ2-7", "TRBJ1-2", "TRBJ2-1"],
                    "c_gene": ["TRBC2", "TRBC1", "TRBC2"],
                    "productive": [True, True, True],
                    "cdr3": ["CASSF", "CASKF", "CASGF"],
                    "cdr3_nt": ["TGTGCCAGC", "TGTGCCAAA", "TGTGCCGGG"],
                    "reads": [50, 30, 10],
                    "umis": [5, 3, 1],
                    "raw_clonotype_id": ["clonotype1", "clonotype2", "clonotype3"],
                }
            ).to_csv(path, index=False)
            result = read_10x(directory)
        self.assertEqual(len(result), 2)
        self.assertEqual(set(result["cell_id"]), {"AA-1", "BB-1"})
        self.assertEqual(result["source"].unique().tolist(), ["10x_contig"])

    def test_10x_clonotypes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "clonotypes.csv"
            pd.DataFrame(
                {
                    "clonotype_id": ["clonotype1"],
                    "frequency": [12],
                    "proportion": [0.4],
                    "cdr3s_aa": ["TRA:CAVR;TRA:CAVS;TRB:CASSF"],
                    "cdr3s_nt": ["TRA:TGTGCT;TRA:TGTGCA;TRB:TGTGCCAGC"],
                }
            ).to_csv(path, index=False)
            result = read_10x(path, sample_id="sc1")
        self.assertEqual(set(result["locus"]), {"TRA", "TRB"})
        self.assertEqual(len(result), 3)
        self.assertEqual(result["source_clonotype_id"].nunique(), 1)

    def test_airr(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "airr_rearrangement.tsv"
            pd.DataFrame(
                {
                    "cell_id": ["AA-1"],
                    "clone_id": ["clonotype1"],
                    "sequence_id": ["AA_contig_1"],
                    "productive": [True],
                    "v_call": ["TRBV7-9*01"],
                    "j_call": ["TRBJ2-7*01"],
                    "junction": ["TGTGCCAGC"],
                    "junction_aa": ["CASSF"],
                    "consensus_count": [50],
                    "duplicate_count": [5],
                    "is_cell": [True],
                }
            ).to_csv(path, sep="\t", index=False)
            result = read_airr(path, sample_id="sc1", cell_only=True)
        self.assertEqual(result.loc[0, "locus"], "TRB")
        self.assertEqual(result.loc[0, "umi_count"], 5)


class AnalysisTests(unittest.TestCase):
    def setUp(self) -> None:
        self.definition = CloneDefinition()
        self.bulk = pd.DataFrame(
            {
                "sample_id": ["pre", "pre", "post", "post"],
                "cell_id": [pd.NA] * 4,
                "sequence_id": ["b1", "b2", "b3", "b4"],
                "source": ["mixcr"] * 4,
                "source_clonotype_id": ["1", "2", "1", "2"],
                "locus": ["TRB"] * 4,
                "productive": [True] * 4,
                "high_confidence": [pd.NA] * 4,
                "junction": ["AAA", "BBB", "AAA", "BBB"],
                "junction_aa": ["K", "N", "K", "N"],
                "v_call": ["TRBV1", "TRBV2", "TRBV1", "TRBV2"],
                "d_call": [pd.NA] * 4,
                "j_call": ["TRBJ1", "TRBJ2", "TRBJ1", "TRBJ2"],
                "c_call": [pd.NA] * 4,
                "read_count": [1, 999, 300, 700],
                "umi_count": [pd.NA] * 4,
                "cell_count": [pd.NA] * 4,
                "frequency": [0.001, 0.999, 0.3, 0.7],
            }
        )
        self.sc = self.bulk.iloc[[2]].copy()
        self.sc["sample_id"] = "sc_post"
        self.sc["cell_id"] = "cell-1"
        self.sc["source"] = "10x_contig"
        self.sc["source_clonotype_id"] = "clonotype1"
        self.sc["read_count"] = 20
        self.sc["umi_count"] = 4

    def test_abundance_link_and_expansion(self) -> None:
        abundance = clone_abundance(self.bulk, self.definition)
        links = link_bulk_to_single_cell(self.bulk, self.sc, self.definition)
        result = test_longitudinal_expansion(
            abundance,
            baseline_sample="pre",
            comparison_samples=["post"],
            fold_change=2,
        )
        matched = links[links["matched"]]
        self.assertEqual(len(matched), 2)  # same clone is present in two bulk samples
        self.assertTrue(matched["clone_id"].str.contains("AAA").all())
        expanded = result[result["clone_id"].str.contains("AAA")].iloc[0]
        self.assertGreater(expanded["fold_change"], 100)
        self.assertTrue(expanded["expanded"])

    def test_phenotype_composition_and_flow(self) -> None:
        sc = pd.concat(
            [
                self.sc.assign(sample_id="t0", cell_id="c1", junction="AAA"),
                self.sc.assign(sample_id="t0", cell_id="c2", junction="AAA"),
                self.sc.assign(sample_id="t1", cell_id="c3", junction="AAA"),
                self.sc.assign(sample_id="t1", cell_id="c4", junction="AAA"),
            ],
            ignore_index=True,
        )
        metadata = pd.DataFrame(
            {
                "sample_id": ["t0", "t0", "t1", "t1"],
                "cell_id": ["c1", "c2", "c3", "c4"],
                "state": ["effector", "effector", "memory", "memory"],
            }
        )
        composition = phenotype_composition(
            sc, metadata, phenotype_col="state", definition=self.definition
        )
        flux = phenotypic_flux(composition, from_sample="t0", to_sample="t1")
        flow = phenotype_flow(composition, from_sample="t0", to_sample="t1")
        self.assertAlmostEqual(float(flux.loc[0, "phenotypic_flux"]), 2.0)
        transition = flow[
            (flow["from_phenotype"] == "effector") & (flow["to_phenotype"] == "memory")
        ]
        self.assertGreater(float(transition["outflow"].sum()), 0)

    def test_project_container(self) -> None:
        project = ImmuneProject(bulk_chains=self.bulk, single_cell_chains=self.sc)
        self.assertFalse(project.build_links().empty)
        self.assertEqual(len(project.repertoire_metrics()), 2)


if __name__ == "__main__":
    unittest.main()
