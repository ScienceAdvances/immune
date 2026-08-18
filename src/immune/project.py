"""High-level container for a longitudinal bulk and single-cell study."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from .bulk import (
    bulk_summary,
    differential_clonotype_abundance,
    gene_usage,
    rarefaction_curve,
    segment_usage,
    skbio_alpha_diversity,
    skbio_beta_diversity,
    spectratype,
)
from .clonotypes import CloneDefinition, clone_abundance, link_bulk_to_single_cell
from .io import read_10x, read_mixcr
from .phenotype import phenotype_composition
from .repertoire import (
    clonality_proportion,
    clonality_rank,
    coverage_diversity,
    hill_diversity,
    public_overlap,
    public_repertoire,
    rank_abundance,
)
from .singlecell import to_scirpy
from .stats import repertoire_metrics, test_longitudinal_expansion


def _empty_chains() -> pd.DataFrame:
    return pd.DataFrame()


@dataclass
class ImmuneProject:
    """A study-level container that keeps bulk and cell-level observations separate."""

    bulk_chains: pd.DataFrame = field(default_factory=_empty_chains)
    single_cell_chains: pd.DataFrame = field(default_factory=_empty_chains)
    sample_metadata: pd.DataFrame = field(default_factory=pd.DataFrame)
    links: pd.DataFrame = field(default_factory=pd.DataFrame)

    def add_mixcr(self, path: str, *, sample_id: str, **kwargs: object) -> ImmuneProject:
        table = read_mixcr(path, sample_id=sample_id, **kwargs)
        self.bulk_chains = pd.concat([self.bulk_chains, table], ignore_index=True)
        return self

    def add_10x(
        self, path: str, *, sample_id: str | None = None, **kwargs: object
    ) -> ImmuneProject:
        table = read_10x(path, sample_id=sample_id, **kwargs)
        self.single_cell_chains = pd.concat(
            [self.single_cell_chains, table], ignore_index=True
        )
        return self

    def bulk_abundance(
        self, definition: CloneDefinition | None = None, *, count_col: str | None = None
    ) -> pd.DataFrame:
        if self.bulk_chains.empty:
            raise ValueError("No bulk repertoire has been added")
        return clone_abundance(self.bulk_chains, definition, count_col=count_col)

    def single_cell_abundance(
        self, definition: CloneDefinition | None = None
    ) -> pd.DataFrame:
        if self.single_cell_chains.empty:
            raise ValueError("No single-cell repertoire has been added")
        return clone_abundance(self.single_cell_chains, definition)

    def build_links(
        self,
        definition: CloneDefinition | None = None,
        *,
        include_unmatched: bool = True,
    ) -> pd.DataFrame:
        if self.bulk_chains.empty or self.single_cell_chains.empty:
            raise ValueError("Both bulk and single-cell receptor data are required")
        self.links = link_bulk_to_single_cell(
            self.bulk_chains,
            self.single_cell_chains,
            definition,
            include_unmatched=include_unmatched,
        )
        return self.links

    def test_expansion(
        self,
        *,
        baseline_sample: str,
        comparison_samples: list[str] | None = None,
        definition: CloneDefinition | None = None,
        **kwargs: object,
    ) -> pd.DataFrame:
        abundance = self.bulk_abundance(definition)
        return test_longitudinal_expansion(
            abundance,
            baseline_sample=baseline_sample,
            comparison_samples=comparison_samples,
            **kwargs,
        )

    def repertoire_metrics(
        self, definition: CloneDefinition | None = None
    ) -> pd.DataFrame:
        return repertoire_metrics(self.bulk_abundance(definition))

    def bulk_diversity(
        self,
        definition: CloneDefinition | None = None,
        **kwargs: object,
    ) -> pd.DataFrame:
        """Compute mature scikit-bio alpha-diversity metrics."""

        return skbio_alpha_diversity(self.bulk_abundance(definition), **kwargs)

    def bulk_distances(
        self,
        definition: CloneDefinition | None = None,
        **kwargs: object,
    ) -> pd.DataFrame:
        """Compute a mature scikit-bio sample distance matrix."""

        return skbio_beta_diversity(self.bulk_abundance(definition), **kwargs)

    def bulk_summary(
        self, definition: CloneDefinition | None = None, **kwargs: object
    ) -> pd.DataFrame:
        """Return Python-native basic and diversity statistics."""

        return bulk_summary(self.bulk_abundance(definition), **kwargs)

    def segment_usage(
        self, definition: CloneDefinition | None = None, **kwargs: object
    ) -> pd.DataFrame:
        return segment_usage(self.bulk_abundance(definition), **kwargs)

    def gene_usage(
        self, definition: CloneDefinition | None = None, **kwargs: object
    ) -> pd.DataFrame:
        return gene_usage(self.bulk_abundance(definition), **kwargs)

    def coverage_diversity(
        self, definition: CloneDefinition | None = None, **kwargs: object
    ) -> pd.DataFrame:
        return coverage_diversity(self.bulk_abundance(definition), **kwargs)

    def hill_diversity(
        self, definition: CloneDefinition | None = None, **kwargs: object
    ) -> pd.DataFrame:
        return hill_diversity(self.bulk_abundance(definition), **kwargs)

    def clonality_proportion(
        self, definition: CloneDefinition | None = None, **kwargs: object
    ) -> pd.DataFrame:
        return clonality_proportion(self.bulk_abundance(definition), **kwargs)

    def clonality_rank(
        self, definition: CloneDefinition | None = None, **kwargs: object
    ) -> pd.DataFrame:
        return clonality_rank(self.bulk_abundance(definition), **kwargs)

    def rank_abundance(
        self, definition: CloneDefinition | None = None, **kwargs: object
    ) -> pd.DataFrame:
        return rank_abundance(self.bulk_abundance(definition), **kwargs)

    def public_repertoire(
        self, definition: CloneDefinition | None = None, **kwargs: object
    ) -> pd.DataFrame:
        return public_repertoire(self.bulk_abundance(definition), **kwargs)

    def public_overlap(
        self, definition: CloneDefinition | None = None, **kwargs: object
    ) -> pd.DataFrame:
        return public_overlap(self.bulk_abundance(definition), **kwargs)

    def spectratype(
        self, definition: CloneDefinition | None = None, **kwargs: object
    ) -> pd.DataFrame:
        return spectratype(self.bulk_abundance(definition), **kwargs)

    def rarefaction(
        self, definition: CloneDefinition | None = None, **kwargs: object
    ) -> pd.DataFrame:
        return rarefaction_curve(self.bulk_abundance(definition), **kwargs)

    def differential_abundance(
        self,
        *,
        design: str,
        contrast: list[str],
        definition: CloneDefinition | None = None,
        **kwargs: object,
    ):
        if self.sample_metadata.empty:
            raise ValueError("sample_metadata is required for differential abundance")
        return differential_clonotype_abundance(
            self.bulk_abundance(definition),
            self.sample_metadata,
            design=design,
            contrast=contrast,
            **kwargs,
        )

    def to_scirpy(self, **kwargs: object):
        """Convert the single-cell receptor table to Scirpy AnnData."""

        if self.single_cell_chains.empty:
            raise ValueError("No single-cell repertoire has been added")
        return to_scirpy(self.single_cell_chains, **kwargs)

    def phenotype_composition(
        self,
        cell_metadata: pd.DataFrame,
        *,
        phenotype_col: str,
        definition: CloneDefinition | None = None,
        **kwargs: object,
    ) -> pd.DataFrame:
        if self.single_cell_chains.empty:
            raise ValueError("No single-cell repertoire has been added")
        return phenotype_composition(
            self.single_cell_chains,
            cell_metadata,
            phenotype_col=phenotype_col,
            definition=definition,
            **kwargs,
        )
