"""Deterministic RNA/VDJ/bulk fixtures with paired donors and ambiguous TRBs."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ._optional import require_dependency
from .schema import canonicalize
from .singlecell import merge_with_transcriptome


def toy_multimodal(*, seed=0, cells_per_sample=24):
    """Return four paired RNA/VDJ samples; RNA also includes VDJ-negative cells."""
    ad = require_dependency("anndata", extra="singlecell", feature="Example RNA data")
    ir = require_dependency("scirpy", extra="singlecell", feature="Example VDJ data")
    from scipy.sparse import csr_matrix

    if cells_per_sample < 8:
        raise ValueError("Use at least eight cells per sample")
    rng = np.random.default_rng(seed)
    rna_obs, airr_cells = [], []
    for sample in range(4):
        for i in range(cells_per_sample):
            cell_id = f"lib{sample}:cell{i}"
            metadata = {
                "sample_id": f"s{sample}",
                "donor_id": f"d{sample // 2}",
                "library_id": f"lib{sample}",
                "barcode": f"cell{i}",
                "timepoint": sample % 2,
                "condition": "pre" if sample % 2 == 0 else "post",
                "cell_type": "T",
                "cell_state": "naive" if (i + sample) % 3 == 0 else "effector",
            }
            rna_obs.append({"cell_id": cell_id, **metadata})
            if i == cells_per_sample - 1:
                continue
            cell = ir.io.AirrCell(cell_id)
            for key in ("sample_id", "donor_id", "library_id", "barcode", "timepoint", "condition"):
                cell[key] = metadata[key]
            clone = i % 3 if sample % 2 == 0 else (0 if i < cells_per_sample // 2 else i % 3)
            sequences = [
                ("TRA", "TGTGCTTTT" if clone == 0 else "TGTGGTTTT", "CAF" if clone == 0 else "CGF"),
                ("TRB", "TGTGCTTTT" if clone < 2 else "TGTACTTTT", "CAF" if clone < 2 else "CTF"),
            ]
            for locus, nt, aa in sequences:
                chain = ir.io.AirrCell.empty_chain_dict()
                chain.update(
                    {
                        "sequence_id": f"{cell_id}_{locus}",
                        "locus": locus,
                        "productive": True,
                        "junction": nt,
                        "junction_aa": aa,
                        "v_call": "TRAV1-1" if locus == "TRA" else "TRBV1",
                        "j_call": "TRAJ1" if locus == "TRA" else "TRBJ1-1",
                        "consensus_count": 20,
                        "duplicate_count": 3,
                    }
                )
                cell.add_chain(chain)
            airr_cells.append(cell)
    obs = pd.DataFrame(rna_obs).set_index("cell_id")
    matrix = rng.poisson(3, size=(len(obs), 40))
    matrix[obs["cell_state"].eq("effector").to_numpy(), 1:5] += 5
    rna = ad.AnnData(
        csr_matrix(matrix),
        obs=obs,
        var=pd.DataFrame(index=["MT-CO1", *[f"G{i}" for i in range(1, 40)]]),
    )
    rna.layers["counts"] = rna.X.copy()
    airr = ir.io.from_airr_cells(airr_cells)
    return merge_with_transcriptome(rna, airr)


def toy_bulk():
    rows = []
    for sample in range(4):
        for junction, aa, count in (
            ("TGTGCTTTT", "CAF", 80 if sample % 2 else 20),
            ("TGTACTTTT", "CTF", 20 if sample % 2 else 80),
        ):
            rows.append(
                {
                    "sample_id": f"s{sample}",
                    "donor_id": f"d{sample // 2}",
                    "locus": "TRB",
                    "junction": junction,
                    "junction_aa": aa,
                    "productive": True,
                    "v_call": "TRBV1",
                    "j_call": "TRBJ1-1",
                    "read_count": count,
                    "source": "example_bulk",
                }
            )
    return canonicalize(pd.DataFrame(rows))
