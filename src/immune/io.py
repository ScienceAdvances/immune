"""Readers for current MiXCR, Cell Ranger V(D)J and AIRR exports."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .schema import (
    SchemaError,
    as_bool,
    canonicalize,
    first_present,
    infer_locus,
    validate_rearrangements,
)


def _read_table(path: str | Path, *, sep: str | None = None) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if sep is None:
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def _series(
    table: pd.DataFrame,
    aliases: list[str],
    *,
    default: object = pd.NA,
) -> pd.Series:
    column = first_present(table.columns, aliases)
    if column is None:
        return pd.Series(default, index=table.index, dtype="object")
    return table[column]


def _sample_ids(table: pd.DataFrame, supplied: str | None, fallback: str) -> pd.Series:
    if supplied is not None:
        return pd.Series(str(supplied), index=table.index, dtype="string")
    values = _series(table, ["sample_id", "sample", "repertoire_id"])
    return values.fillna(fallback).astype("string").replace("", fallback)


def _filter_rows(
    table: pd.DataFrame,
    *,
    productive_only: bool,
    high_confidence_only: bool = False,
    cell_only: bool = False,
) -> pd.DataFrame:
    keep = pd.Series(True, index=table.index)
    if productive_only and table["productive"].notna().any():
        keep &= table["productive"].fillna(False)
    if high_confidence_only and table["high_confidence"].notna().any():
        keep &= table["high_confidence"].fillna(False)
    if cell_only and "is_cell" in table and table["is_cell"].notna().any():
        keep &= table["is_cell"].map(as_bool).fillna(False)
    return table.loc[keep].reset_index(drop=True)


def read_airr(
    path: str | Path,
    *,
    sample_id: str | None = None,
    source: str = "airr",
    productive_only: bool = True,
    cell_only: bool = False,
) -> pd.DataFrame:
    """Read an AIRR rearrangement TSV into the immune canonical schema."""

    path = Path(path)
    raw = _read_table(path, sep="\t")
    required_any = {"junction", "junction_aa", "sequence_id", "v_call", "j_call"}
    if not required_any.intersection(raw.columns):
        raise SchemaError(f"{path} does not look like an AIRR rearrangement table")

    result = pd.DataFrame(index=raw.index)
    result["sample_id"] = _sample_ids(raw, sample_id, path.stem)
    result["cell_id"] = _series(raw, ["cell_id", "cell", "barcode"])
    result["sequence_id"] = _series(raw, ["sequence_id", "contig_id"])
    result["source"] = source
    result["source_clonotype_id"] = _series(raw, ["clone_id", "clonotype_id"])
    result["locus"] = _series(raw, ["locus", "chain"])
    result["productive"] = _series(raw, ["productive"])
    result["high_confidence"] = _series(raw, ["high_confidence"])
    result["junction"] = _series(raw, ["junction", "cdr3_nt"])
    result["junction_aa"] = _series(raw, ["junction_aa", "cdr3"])
    result["v_call"] = _series(raw, ["v_call", "v_gene"])
    result["d_call"] = _series(raw, ["d_call", "d_gene"])
    result["j_call"] = _series(raw, ["j_call", "j_gene"])
    result["c_call"] = _series(raw, ["c_call", "c_gene"])
    result["read_count"] = _series(raw, ["consensus_count", "read_count", "reads"])
    result["umi_count"] = _series(raw, ["duplicate_count", "umi_count", "umis"])
    result["cell_count"] = _series(raw, ["cell_count"])
    result["frequency"] = _series(raw, ["frequency", "clone_fraction"])
    if "is_cell" in raw:
        result["is_cell"] = raw["is_cell"]

    result = canonicalize(result)
    result = _filter_rows(result, productive_only=productive_only, cell_only=cell_only)
    validate_rearrangements(result)
    return result


def _read_10x_contigs(
    path: Path,
    *,
    sample_id: str | None,
    productive_only: bool,
    high_confidence_only: bool,
) -> pd.DataFrame:
    raw = _read_table(path, sep=",")
    if "barcode" not in raw or "cdr3_nt" not in raw:
        raise SchemaError(f"{path} is not a Cell Ranger contig annotations CSV")

    result = pd.DataFrame(index=raw.index)
    result["sample_id"] = _sample_ids(raw, sample_id, path.parent.name)
    result["cell_id"] = raw["barcode"]
    result["sequence_id"] = _series(raw, ["contig_id"])
    result["source"] = "10x_contig"
    result["source_clonotype_id"] = _series(
        raw, ["raw_clonotype_id", "clonotype_id", "exact_subclonotype_id"]
    )
    result["locus"] = _series(raw, ["chain", "locus"])
    result["productive"] = _series(raw, ["productive"])
    result["high_confidence"] = _series(raw, ["high_confidence"])
    result["junction"] = raw["cdr3_nt"]
    result["junction_aa"] = _series(raw, ["cdr3"])
    result["v_call"] = _series(raw, ["v_gene"])
    result["d_call"] = _series(raw, ["d_gene"])
    result["j_call"] = _series(raw, ["j_gene"])
    result["c_call"] = _series(raw, ["c_gene"])
    result["read_count"] = _series(raw, ["reads"])
    result["umi_count"] = _series(raw, ["umis"])
    result["cell_count"] = pd.NA
    result["frequency"] = pd.NA
    if "is_cell" in raw:
        result["is_cell"] = raw["is_cell"]

    result = canonicalize(result)
    result = _filter_rows(
        result,
        productive_only=productive_only,
        high_confidence_only=high_confidence_only,
        cell_only=True,
    )
    validate_rearrangements(result)
    return result


def _parse_chain_pairs(value: object) -> list[tuple[str, str]]:
    if value is None or pd.isna(value):
        return []
    pairs: list[tuple[str, str]] = []
    for item in str(value).split(";"):
        if ":" not in item:
            continue
        chain, sequence = item.split(":", 1)
        if chain and sequence:
            pairs.append((chain.upper(), sequence.upper()))
    return pairs


def _read_10x_clonotypes(path: Path, *, sample_id: str | None) -> pd.DataFrame:
    raw = _read_table(path, sep=",")
    required = {"clonotype_id", "frequency", "cdr3s_nt"}
    if not required.issubset(raw.columns):
        raise SchemaError(f"{path} is not a current Cell Ranger clonotypes.csv file")

    records: list[dict[str, object]] = []
    fallback_sample = path.parent.name
    for _, row in raw.iterrows():
        nt_pairs = _parse_chain_pairs(row.get("cdr3s_nt"))
        aa_pairs = _parse_chain_pairs(row.get("cdr3s_aa"))
        chains = sorted({chain for chain, _ in nt_pairs} | {chain for chain, _ in aa_pairs})
        for chain in chains:
            nt_sequences = [sequence for pair_chain, sequence in nt_pairs if pair_chain == chain]
            aa_sequences = [sequence for pair_chain, sequence in aa_pairs if pair_chain == chain]
            for occurrence in range(max(len(nt_sequences), len(aa_sequences))):
                records.append(
                    {
                        "sample_id": sample_id or row.get("sample") or fallback_sample,
                        "cell_id": pd.NA,
                        "sequence_id": f"{row['clonotype_id']}:{chain}:{occurrence + 1}",
                        "source": "10x_clonotype",
                        "source_clonotype_id": row["clonotype_id"],
                        "locus": chain,
                        "productive": True,
                        "high_confidence": True,
                        "junction": nt_sequences[occurrence]
                        if occurrence < len(nt_sequences)
                        else pd.NA,
                        "junction_aa": aa_sequences[occurrence]
                        if occurrence < len(aa_sequences)
                        else pd.NA,
                        "read_count": pd.NA,
                        "umi_count": pd.NA,
                        "cell_count": row.get("frequency", pd.NA),
                        "frequency": row.get("proportion", pd.NA),
                    }
                )
    result = canonicalize(pd.DataFrame.from_records(records))
    validate_rearrangements(result)
    return result


def read_10x(
    path: str | Path,
    *,
    sample_id: str | None = None,
    productive_only: bool = True,
    high_confidence_only: bool = True,
) -> pd.DataFrame:
    """Read current Cell Ranger V(D)J output.

    ``path`` may point to ``airr_rearrangement.tsv``,
    ``filtered_contig_annotations.csv``, ``all_contig_annotations.csv``,
    ``clonotypes.csv``, or a directory containing one of them. Directories prefer
    cell-level filtered contigs, then AIRR, then the high-level clonotype table.
    """

    path = Path(path)
    if path.is_dir():
        candidates = [
            path / "filtered_contig_annotations.csv",
            path / "airr_rearrangement.tsv",
            path / "all_contig_annotations.csv",
            path / "clonotypes.csv",
        ]
        try:
            path = next(candidate for candidate in candidates if candidate.exists())
        except StopIteration as exc:
            raise FileNotFoundError(f"No supported Cell Ranger V(D)J output found in {path}") from exc

    if path.name == "airr_rearrangement.tsv" or path.suffix.lower() == ".tsv":
        return read_airr(
            path,
            sample_id=sample_id,
            source="10x_airr",
            productive_only=productive_only,
            cell_only=True,
        )
    if path.name == "clonotypes.csv":
        return _read_10x_clonotypes(path, sample_id=sample_id)
    return _read_10x_contigs(
        path,
        sample_id=sample_id,
        productive_only=productive_only,
        high_confidence_only=high_confidence_only,
    )


def read_mixcr(
    path: str | Path,
    *,
    sample_id: str | None = None,
    productive_only: bool = True,
) -> pd.DataFrame:
    """Read a MiXCR ``exportClones`` table or MiXCR AIRR export.

    MiXCR export columns depend on the preset and explicit export flags. This
    reader recognizes the common default fields and the explicit ``-count``,
    ``-readCount``, ``-uniqueTagCount Molecule``, ``-vHit``, ``-jHit``,
    ``-nFeature CDR3`` and ``-aaFeature CDR3`` forms.
    """

    path = Path(path)
    raw = _read_table(path)
    if {"junction", "v_call", "j_call"}.intersection(raw.columns):
        return read_airr(
            path,
            sample_id=sample_id,
            source="mixcr_airr",
            productive_only=productive_only,
        )

    nt_column = first_present(
        raw.columns,
        ["nSeqCDR3", "nSeqImputedCDR3", "nFeatureCDR3", "cdr3_nt", "junction"],
    )
    aa_column = first_present(
        raw.columns,
        ["aaSeqCDR3", "aaSeqImputedCDR3", "aaFeatureCDR3", "cdr3_aa", "junction_aa"],
    )
    if nt_column is None and aa_column is None:
        raise SchemaError(
            "MiXCR table has no recognized CDR3 field; export nSeqCDR3 and/or aaSeqCDR3"
        )

    result = pd.DataFrame(index=raw.index)
    result["sample_id"] = _sample_ids(raw, sample_id, path.stem)
    result["cell_id"] = _series(raw, ["tagValueCELL", "cellId", "cell_id"])
    result["sequence_id"] = _series(raw, ["cloneId", "clone_id", "sequence_id"])
    result["source"] = "mixcr"
    result["source_clonotype_id"] = _series(raw, ["cloneId", "clone_id", "cellGroup"])
    result["locus"] = _series(raw, ["chains", "topChains", "chain", "locus"])
    result["productive"] = _series(raw, ["isProductive", "productive"])
    result["high_confidence"] = pd.NA
    result["junction"] = raw[nt_column] if nt_column is not None else pd.NA
    result["junction_aa"] = raw[aa_column] if aa_column is not None else pd.NA
    result["v_call"] = _series(
        raw, ["bestVHit", "allVHitsWithScore", "vHit", "vGene", "bestVGene", "v_call"]
    )
    result["d_call"] = _series(
        raw, ["bestDHit", "allDHitsWithScore", "dHit", "dGene", "bestDGene", "d_call"]
    )
    result["j_call"] = _series(
        raw, ["bestJHit", "allJHitsWithScore", "jHit", "jGene", "bestJGene", "j_call"]
    )
    result["c_call"] = _series(
        raw, ["bestCHit", "allCHitsWithScore", "cHit", "cGene", "bestCGene", "c_call"]
    )
    result["read_count"] = _series(raw, ["cloneCount", "readCount", "count"])
    result["umi_count"] = _series(raw, ["uniqueTagCountMolecule", "moleculeCount"])
    result["cell_count"] = _series(raw, ["uniqueTagCountCell", "cellCount"])
    result["frequency"] = _series(raw, ["cloneFraction", "readFraction", "fraction"])

    if result["productive"].isna().all():
        aa = result["junction_aa"].astype("string")
        result["productive"] = aa.map(
            lambda value: pd.NA
            if pd.isna(value)
            else not bool(pd.Series([value]).str.contains(r"\*|_", regex=True).iloc[0])
        )

    result = canonicalize(result)
    for index, row in result[result["locus"].isna()].iterrows():
        result.at[index, "locus"] = infer_locus(
            row.get("v_call"), row.get("d_call"), row.get("j_call"), row.get("c_call")
        )
    result = _filter_rows(result, productive_only=productive_only)
    validate_rearrangements(result)
    return result
