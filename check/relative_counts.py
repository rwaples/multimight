"""Check a converter output against pedigree-graph, person by person.

For every proband, trait and kind of relative in the pool, relatives_n and
relatives_n_trait must equal the counts pedigree-graph's classifier gives on
the same pedigree. Exits 1 on any difference.

Usage:
    pixi run check REP_DIR OUT.tsv.gz [NAME=N ...]

REP_DIR holds the replicate's cohort.parquet and pedigree.parquet. A trait
named traitN in the pool is simACE trait N; for other names (from the
converter's --trait-names) pass NAME=N, e.g. early=1 late=2.
"""

from __future__ import annotations

import sys

import numpy as np
import polars as pl
from pedigree_graph import PedigreeGraph

# The converter's kinds, as pedigree-graph registry codes
KINDS = {
    "PO": ("MO", "FO"),
    "FS": ("FS",),
    "HS": ("MHS", "PHS"),
    "mHS": ("MHS",),
    "pHS": ("PHS",),
    "Av": ("Av",),
    "1G": ("GP",),
    "1C": ("1C",),
}


def main(rep: str, out: str, names: dict[str, int]) -> int:
    """Compare the pool at ``out`` with pedigree-graph's counts for ``rep``; return the exit code."""
    pool = pl.read_csv(out, separator="\t", schema_overrides={"person_id": pl.Int64})
    kinds = pool["relatives_kind"].unique(maintain_order=True).to_list()
    unknown = sorted(set(kinds) - set(KINDS))
    if unknown:
        print(f"unknown relatives_kind in {out}: {unknown}")
        return 1
    traits = pool["trait"].unique(maintain_order=True).to_list()
    for trait in traits:
        if trait not in names and trait.startswith("trait") and trait[5:].isdigit():
            names[trait] = int(trait[5:])
    missing = [t for t in traits if t not in names]
    if missing:
        print(f"which simACE trait is {missing}? pass NAME=N")
        return 1

    ped = pl.read_parquet(f"{rep}/pedigree.parquet", columns=["id", "generation", "sex", "mother", "father", "twin"])
    trait_cols = sorted({f"affected{n}" for n in names.values()})
    cohort = pl.read_parquet(f"{rep}/cohort.parquet", columns=["id", *trait_cols])
    cohort = cohort.filter(pl.col(trait_cols[0]).is_not_null())
    df = ped.join(cohort, on="id", how="left")
    phenotyped = df[trait_cols[0]].is_not_null().to_numpy()

    thresholds = {"phenotyped": (np.where(phenotyped, 0.0, np.nan), 0.0)}
    for col in trait_cols:
        thresholds[col] = (np.where(df[col].fill_null(False).to_numpy(), 0.0, np.nan), 0.0)
    codes = tuple(dict.fromkeys(code for kind in kinds for code in KINDS[kind]))
    counts = PedigreeGraph.from_frame(df).relatives_per_person(categories=codes, thresholds=thresholds, progress=False)
    ids = df["id"].to_numpy().astype(np.int64)

    failed = False
    for trait in traits:
        rows = pool.filter(pl.col("trait") == trait)
        affected_col = f"affected{names[trait]}"
        for kind in kinds:
            expected = pl.DataFrame(
                {
                    "person_id": ids,
                    "n_expected": counts.sum(KINDS[kind], column="phenotyped").astype(np.int64),
                    "n_trait_expected": counts.sum(KINDS[kind], column=affected_col).astype(np.int64),
                }
            )
            block = rows.filter(pl.col("relatives_kind") == kind).join(expected, on="person_id", how="left")
            bad = block.filter(
                (pl.col("relatives_n") != pl.col("n_expected"))
                | (pl.col("relatives_n_trait") != pl.col("n_trait_expected"))
                | pl.col("n_expected").is_null()
            )
            failed |= bad.height > 0
            print(f"{trait:10s} {kind:4s} {block.height:9d} probands  {bad.height:7d} differ")
            if bad.height:
                print(bad.head(5))
    print("FAIL" if failed else "PASS: every count matches pedigree-graph")
    return 1 if failed else 0


if __name__ == "__main__":
    if len(sys.argv) < 3 or not all("=" in a for a in sys.argv[3:]):
        sys.exit(__doc__)
    names = {name: int(n) for name, n in (a.split("=", 1) for a in sys.argv[3:])}
    sys.exit(main(sys.argv[1], sys.argv[2], names))
