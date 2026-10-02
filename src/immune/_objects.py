"""Native scverse object contract; no dependency on cellscope or Muon."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

import pandas as pd

from ._optional import require_dependency


def airr_data(data, airr_mod="airr"):
    adata = data.mod[airr_mod] if hasattr(data, "mod") else data
    if not hasattr(adata, "obs") or not hasattr(adata, "obsm"):
        raise TypeError("Expected AIRR AnnData or MuData")
    if not adata.obs_names.is_unique:
        raise ValueError("Cell identifiers must be unique")
    if "airr" not in adata.obsm:
        raise ValueError("AIRR records are required in obsm['airr']")
    return adata


def record(data, key, table=None, *, params=None, backend="immune"):
    try:
        backend_version = version(backend)
    except PackageNotFoundError:
        backend_version = "source"
    entry = {
        "params": params or {},
        "backend": backend,
        "backend_version": backend_version,
        "schema_version": "1.0",
    }
    if table is not None:
        entry["table"] = table
    data.uns.setdefault("immune", {})[key] = entry


def sync_obs(data, airr_mod="airr"):
    if hasattr(data, "mod"):
        for col in data.mod[airr_mod].obs:
            data.obs[f"{airr_mod}:{col}"] = data.mod[airr_mod].obs[col].reindex(data.obs_names)


def cell_obs(data, *, airr_mod="airr", gex_mod="gex", clone_key="clone_id"):
    """Resolve receptor annotations and aligned RNA labels without matrix copies."""
    adata = airr_data(data, airr_mod)
    if clone_key not in adata.obs:
        raise KeyError(f"Run define_clonotypes before accessing {clone_key!r}")
    frame = adata.obs.copy()
    if hasattr(data, "mod"):
        if gex_mod not in data.mod:
            raise KeyError(gex_mod)
        gex = data.mod[gex_mod]
        if not gex.obs_names.is_unique:
            raise ValueError("RNA cell identifiers must be unique")
        for col in gex.obs:
            values = gex.obs[col].reindex(frame.index)
            if col not in frame:
                frame[col] = values
            elif col in {"sample_id", "donor_id", "library_id"}:
                conflict = frame[col].notna() & values.notna() & frame[col].ne(values)
                if conflict.any():
                    raise ValueError(f"Conflicting {col} between modalities")
                frame[col] = frame[col].combine_first(values)
            frame[f"{gex_mod}:{col}"] = values
        frame["has_gex"] = frame.index.isin(gex.obs_names)
    frame["cell_id"] = frame.index.astype(str)
    return frame


def merge_modalities(rna, airr, *, rna_mod="gex", airr_mod="airr", join="outer"):
    md = require_dependency("mudata", extra="singlecell", feature="RNA/VDJ integration")
    if rna_mod == airr_mod or join not in {"inner", "outer"}:
        raise ValueError("Use distinct modality names and join='inner' or 'outer'")
    for obj in (rna, airr):
        if not obj.obs_names.is_unique:
            raise ValueError("Cell IDs must be unique; assign library:barcode IDs before merging")
    for key in ("sample_id", "donor_id", "library_id", "barcode", "condition", "timepoint"):
        if key in rna.obs and key in airr.obs:
            compared = pd.concat([rna.obs[key], airr.obs[key]], axis=1)
            if compared.nunique(axis=1, dropna=True).gt(1).any():
                raise ValueError(f"Conflicting {key} for shared cell IDs")
    if join == "inner":
        common = rna.obs_names.intersection(airr.obs_names, sort=False)
        rna, airr = rna[common].copy(), airr[common].copy()
    data = md.MuData({rna_mod: rna, airr_mod: airr})
    data.pull_obs()
    for mod, obj in ((rna_mod, rna), (airr_mod, airr)):
        for column in obj.obs:
            data.obs[f"{mod}:{column}"] = obj.obs[column].reindex(data.obs_names)
    data.obs["has_gex"] = data.obs_names.isin(rna.obs_names)
    data.obs["has_airr"] = data.obs_names.isin(airr.obs_names)
    for key in ("sample_id", "donor_id", "library_id", "barcode", "condition", "timepoint"):
        columns = [obj.obs[key].reindex(data.obs_names) for obj in (rna, airr) if key in obj.obs]
        if columns:
            data.obs[key] = pd.concat(columns, axis=1).bfill(axis=1).iloc[:, 0]
    record(
        data, "merge_modalities", params={"join": join, "gex_mod": rna_mod, "airr_mod": airr_mod}
    )
    return data
