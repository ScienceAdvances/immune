"""Clonotype definitions, abundance tables and cross-platform linking."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .schema import clean_gene_call, validate_rearrangements


@dataclass(frozen=True, slots=True)
class CloneDefinition:
    """Immutable rules for assigning clonotype identifiers.

    The default is deliberately conservative for bulk-to-single-cell matching:
    nucleotide CDR3 on the beta chain without V/J calls. V/J calls can differ
    between annotation pipelines even when the rearrangement is identical.
    """

    sequence: str = "junction"
    loci: tuple[str, ...] = ("TRB",)
    use_v: bool = False
    use_j: bool = False
    use_alleles: bool = False

    def __post_init__(self) -> None:
        if self.sequence not in {"junction", "junction_aa"}:
            raise ValueError("sequence must be 'junction' or 'junction_aa'")
        if not self.loci:
            raise ValueError("At least one locus is required")

    @property
    def label(self) -> str:
        sequence = "NT" if self.sequence == "junction" else "AA"
        genes = "".join(["+V" if self.use_v else "", "+J" if self.use_j else ""])
        return f"{sequence}{genes}:{','.join(self.loci)}"


def _normalize_gene_for_key(value: object, use_alleles: bool) -> str:
    gene = clean_gene_call(value, keep_allele=use_alleles)
    return "?" if pd.isna(gene) else str(gene)


def _clone_key(row: pd.Series, definition: CloneDefinition) -> object:
    locus = str(row.get("locus", "")).upper()
    if locus not in definition.loci:
        return pd.NA
    sequence = row.get(definition.sequence)
    if sequence is None or pd.isna(sequence) or str(sequence).strip() == "":
        return pd.NA
    parts = [locus, "NT" if definition.sequence == "junction" else "AA", str(sequence).upper()]
    if definition.use_v:
        v_call = _normalize_gene_for_key(row.get("v_call"), definition.use_alleles)
        if v_call == "?":
            return pd.NA
        parts.extend(["V", v_call])
    if definition.use_j:
        j_call = _normalize_gene_for_key(row.get("j_call"), definition.use_alleles)
        if j_call == "?":
            return pd.NA
        parts.extend(["J", j_call])
    return "|".join(parts)


def add_clonotype_ids(
    chains: pd.DataFrame,
    definition: CloneDefinition | None = None,
    *,
    key_added: str = "clone_id",
) -> pd.DataFrame:
    """Assign stable clonotype keys to chain-level rearrangements."""

    definition = definition or CloneDefinition()
    validate_rearrangements(chains)
    result = chains.copy()
    result[key_added] = result.apply(lambda row: _clone_key(row, definition), axis=1).astype("string")
    result["clone_definition"] = definition.label
    return result


def _choose_count_column(table: pd.DataFrame, requested: str | None) -> tuple[str, str]:
    if requested is not None:
        if requested not in table.columns:
            raise KeyError(f"Count column {requested!r} is not present")
        return requested, requested.removesuffix("_count")
    for column, unit in (
        ("cell_count", "cell"),
        ("umi_count", "umi"),
        ("read_count", "read"),
    ):
        if column in table and table[column].notna().any():
            return column, unit
    raise ValueError("No cell_count, umi_count or read_count is available")


def clone_abundance(
    chains: pd.DataFrame,
    definition: CloneDefinition | None = None,
    *,
    count_col: str | None = None,
) -> pd.DataFrame:
    """Collapse chain-level data to sample-by-clonotype abundance.

    Cell-level tables are counted by unique ``cell_id``. Bulk and high-level
    clonotype tables use an explicit count column, or choose cell, UMI and read
    counts in that order.
    """

    definition = definition or CloneDefinition()
    assigned = add_clonotype_ids(chains, definition)
    assigned = assigned[assigned["clone_id"].notna()].copy()
    if assigned.empty:
        return pd.DataFrame(
            columns=[
                "sample_id",
                "clone_id",
                "locus",
                "junction",
                "junction_aa",
                "v_call",
                "d_call",
                "j_call",
                "c_call",
                "count",
                "count_unit",
                "frequency",
            ]
        )

    metadata_columns = [
        "locus",
        "junction",
        "junction_aa",
        "v_call",
        "d_call",
        "j_call",
        "c_call",
    ]
    sample_results: list[pd.DataFrame] = []
    for sample_id, sample in assigned.groupby("sample_id", observed=True, sort=False):
        sample_keys = ["clone_id"]
        if sample["cell_id"].notna().any():
            grouped_count = (
                sample.dropna(subset=["cell_id"])
                .groupby(sample_keys, observed=True)["cell_id"]
                .nunique()
                .rename("count")
            )
            unit = "cell"
        else:
            selected, unit = _choose_count_column(sample, count_col)
            # A high-level paired clonotype has one row per chain. The selected
            # definition normally filters to one locus; max avoids double
            # counting duplicated rows carrying the same source abundance.
            grouped_count = (
                sample.groupby(sample_keys, observed=True)[selected]
                .max()
                .fillna(0)
                .rename("count")
            )
        metadata = sample.groupby(sample_keys, observed=True)[metadata_columns].first()
        sample_result = metadata.join(grouped_count).reset_index()
        sample_result.insert(0, "sample_id", sample_id)
        sample_result["count_unit"] = unit
        sample_results.append(sample_result)

    result = pd.concat(sample_results, ignore_index=True)
    result["count"] = pd.to_numeric(result["count"], errors="coerce").fillna(0.0)
    totals = result.groupby("sample_id", observed=True)["count"].transform("sum")
    result["frequency"] = np.divide(
        result["count"],
        totals,
        out=np.zeros(len(result), dtype=float),
        where=totals.to_numpy() > 0,
    )
    return result.sort_values(["sample_id", "count"], ascending=[True, False]).reset_index(drop=True)


def _source_ids(values: pd.Series) -> str:
    unique = sorted({str(value) for value in values.dropna() if str(value)})
    return ";".join(unique)


def link_bulk_to_single_cell(
    bulk_chains: pd.DataFrame,
    single_cell_chains: pd.DataFrame,
    definition: CloneDefinition | None = None,
    *,
    include_unmatched: bool = True,
) -> pd.DataFrame:
    """Link bulk clonotypes to single-cell clonotypes under explicit rules."""

    definition = definition or CloneDefinition()
    bulk = add_clonotype_ids(bulk_chains, definition)
    cells = add_clonotype_ids(single_cell_chains, definition)
    bulk = bulk[bulk["clone_id"].notna()].copy()
    cells = cells[cells["clone_id"].notna()].copy()

    bulk_counts = clone_abundance(bulk, definition).rename(
        columns={
            "sample_id": "bulk_sample_id",
            "count": "bulk_count",
            "frequency": "bulk_frequency",
            "count_unit": "bulk_count_unit",
        }
    )
    bulk_sources = (
        bulk.groupby(["sample_id", "clone_id"], observed=True)["source_clonotype_id"]
        .agg(_source_ids)
        .rename("bulk_source_clonotype_ids")
        .reset_index()
        .rename(columns={"sample_id": "bulk_sample_id"})
    )
    bulk_counts = bulk_counts.merge(bulk_sources, on=["bulk_sample_id", "clone_id"], how="left")

    sc_counts = clone_abundance(cells, definition)[
        ["sample_id", "clone_id", "count", "count_unit"]
    ].rename(columns={"count": "n_sc_cells", "count_unit": "sc_count_unit"})
    sc_sources = (
        cells.groupby(["sample_id", "clone_id"], observed=True)["source_clonotype_id"]
        .agg(_source_ids)
        .rename("sc_source_clonotype_ids")
        .reset_index()
    )
    sc = sc_counts.merge(sc_sources, on=["sample_id", "clone_id"], how="left").rename(
        columns={"sample_id": "sc_sample_id"}
    )
    sc["n_sc_source_clonotypes"] = sc["sc_source_clonotype_ids"].fillna("").map(
        lambda value: len([item for item in value.split(";") if item])
    )

    how = "left" if include_unmatched else "inner"
    links = bulk_counts.merge(sc, on="clone_id", how=how)
    links["matched"] = links["sc_sample_id"].notna()
    links["match_type"] = np.where(
        links["matched"],
        "exact_nt_vj"
        if definition.sequence == "junction" and (definition.use_v or definition.use_j)
        else "exact_nt"
        if definition.sequence == "junction"
        else "exact_aa",
        "unmatched",
    )
    links["ambiguous"] = links["n_sc_source_clonotypes"].fillna(0).gt(1)
    links["clone_definition"] = definition.label
    return links.reset_index(drop=True)
