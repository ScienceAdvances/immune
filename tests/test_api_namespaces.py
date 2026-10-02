from __future__ import annotations

import unittest
from unittest.mock import patch

import immune as iu


class NamespaceTests(unittest.TestCase):
    def test_scanpy_style_namespaces_are_public(self) -> None:
        self.assertIs(iu.io.read_mixcr, iu.read_mixcr)
        self.assertIs(iu.pp.clone_abundance, iu.clone_abundance)
        self.assertIs(iu.tl.bulk_summary, iu.bulk_summary)
        self.assertIs(iu.tl.alpha_diversity, iu.skbio_alpha_diversity)
        self.assertIs(iu.pl.bulk_diversity, iu.plot_bulk_diversity)
        self.assertEqual(iu.__version__, "0.5.0")

    def test_single_cell_plot_wrappers_dispatch_to_scirpy(self) -> None:
        data = object()
        with patch("immune.pl.plot_scirpy", return_value="axis") as plot:
            result = iu.pl.clonotype_network(data, color="sample_id")
        self.assertEqual(result, "axis")
        plot.assert_called_once_with(data, "clonotype_network", color="sample_id")


if __name__ == "__main__":
    unittest.main()
