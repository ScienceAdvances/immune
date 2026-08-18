"""Canonical AIRR-like schema and normalization helpers."""

from __future__ import annotations

import re
from collections.abc import Iterable

import numpy as np
import pandas as pd

CANONICAL_COLUMNS = [
    "sample_id",
    "cell_id",
    "sequence_id",
    "source",
    "source_clonotype_id",
    "locus",
    "productive",
    "high_confidence",
    "junction",
    "junction_aa",
    "v_call",
    "d_call",
    "j_call",
    "c_call",
    "read_count",
    "umi_count",
    "cell_count",
    "frequency",
]

NUMERIC_COLUMNS = ["read_count", "umi_count", "cell_count", "frequency"]
BOOLEAN_COLUMNS = ["productive", "high_confidence"]


class SchemaError(ValueError):
    """Raised when a receptor table cannot be converted to the canonical schema."""


def first_present(columns: Iterable[str], aliases: Iterable[str]) -> str | None:
    """Return the first alias present in ``columns`` (case sensitive first)."""

    column_list = list(columns)
    for alias in aliases:
        if alias in column_list:
            return alias
    lower = {str(column).lower(): column for column in column_list}
    for alias in aliases:
        if alias.lower() in lower:
            return lower[alias.lower()]
    return None


def as_bool(value: object) -> object:
    """Parse common AIRR, MiXCR and Cell Ranger boolean encodings."""

    if value is None or pd.isna(value):
        return pd.NA
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float)):
        return bool(value)
    normalized = str(value).strip().lower()
    if normalized in {"true", "t", "1", "yes", "y"}:
        return True
    if normalized in {"false", "f", "0", "no", "n"}:
        return False
    return pd.NA


def clean_sequence(value: object) -> object:
    """Normalize nucleotide/amino-acid sequences while preserving missing values."""

    if value is None or pd.isna(value):
        return pd.NA
    sequence = re.sub(r"\s+", "", str(value)).upper()
    return sequence if sequence else pd.NA


def clean_gene_call(value: object, *, keep_allele: bool = True) -> object:
    """Return the best gene call and remove MiXCR scores and surrounding whitespace."""

    if value is None or pd.isna(value):
        return pd.NA
    gene = str(value).strip()
    if not gene:
        return pd.NA
    gene = re.split(r"[,;]", gene, maxsplit=1)[0]
    gene = re.sub(r"\([^)]*\)$", "", gene).strip()
    if not keep_allele:
        gene = gene.split("*")[0]
    return gene.upper() if gene else pd.NA


def infer_locus(*values: object) -> object:
    """Infer receptor locus from chain or V/D/J/C gene calls."""

    prefixes = ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL")
    for value in values:
        if value is None or pd.isna(value):
            continue
        text = str(value).upper()
        for prefix in prefixes:
            if prefix in text:
                return prefix
    return pd.NA


def canonicalize(table: pd.DataFrame) -> pd.DataFrame:
    """Add missing canonical columns and normalize common data types."""

    result = table.copy()
    for column in CANONICAL_COLUMNS:
        if column not in result:
            result[column] = pd.NA

    for column in ("junction", "junction_aa"):
        result[column] = result[column].map(clean_sequence)
    for column in ("v_call", "d_call", "j_call", "c_call"):
        result[column] = result[column].map(clean_gene_call)
    for column in BOOLEAN_COLUMNS:
        result[column] = result[column].map(as_bool).astype("boolean")
    for column in NUMERIC_COLUMNS:
        result[column] = pd.to_numeric(result[column], errors="coerce")

    missing_locus = result["locus"].isna() | result["locus"].astype("string").str.strip().eq("")
    if missing_locus.any():
        inferred = result.loc[missing_locus].apply(
            lambda row: infer_locus(
                row.get("locus"),
                row.get("v_call"),
                row.get("d_call"),
                row.get("j_call"),
                row.get("c_call"),
            ),
            axis=1,
        )
        result.loc[missing_locus, "locus"] = inferred
    result["locus"] = result["locus"].astype("string").str.upper()

    result["sample_id"] = result["sample_id"].astype("string")
    return result[CANONICAL_COLUMNS + [c for c in result.columns if c not in CANONICAL_COLUMNS]]


def validate_rearrangements(
    table: pd.DataFrame,
    *,
    require_junction: bool = True,
    require_sample: bool = True,
) -> None:
    """Validate the minimum information required by downstream analyses."""

    missing = [column for column in CANONICAL_COLUMNS if column not in table.columns]
    if missing:
        raise SchemaError(f"Table is not canonical; missing columns: {missing}")
    if require_sample and (table["sample_id"].isna().all() or table["sample_id"].eq("").all()):
        raise SchemaError("At least one non-empty sample_id is required")
    if require_junction and table["junction"].isna().all() and table["junction_aa"].isna().all():
        raise SchemaError("No CDR3/junction nucleotide or amino-acid sequence was found")
