"""Compare pipeline outputs against a golden run.

A manifest TSV lists `golden_path  new_path  mode  [key_regex]`, relative to the
two run directories. `{prefix}` in golden paths is replaced by --prefix. A `{}`
placeholder expands to every golden match, so one line covers a directory of
per-sample files; key_regex, if given, maps the golden key to the new one
(first capture group). Missing new files are counted, not failed, so partially
built runs can be checked.

Modes:
  exact   byte-identical; if not, reports whether only row order differs.
  rows    same lines in any order.
  sets    as exact, after sorting the items of every ';'-joined column.
  zscore  headerless anchor/zscore/... table; same anchors, |dz| <= --tol.
  maxz    anchor/max_zscore table with header; same anchors, |dz| <= --tol.
"""

import argparse
import filecmp
import re
import subprocess
import sys
from pathlib import Path

import polars as pl


def build_parser():
    """Return the argument parser for compare.py."""
    p = argparse.ArgumentParser(description="Compare outputs against a golden run.")
    p.add_argument("golden", type=Path, help="Golden run directory.")
    p.add_argument("new", type=Path, help="New run directory.")
    p.add_argument("-m", "--manifest", type=Path, required=True, help="Manifest TSV.")
    p.add_argument("--prefix", default="CUSTOM", help="Golden file prefix (CUSTOM, ENCODE, MERGED).")
    p.add_argument("--only", default=None, help="Regex; check matching manifest lines only.")
    p.add_argument("--tol", type=float, default=1e-12, help="Absolute z-score tolerance.")
    return p


def _sorted_equal(a, b):
    """True if the two files hold the same lines in any order."""
    cmd = "cmp -s <(LC_ALL=C sort {0}) <(LC_ALL=C sort {1})"
    return subprocess.run(
        ["bash", "-c", cmd.format(f"'{a}'", f"'{b}'")]
    ).returncode == 0


def _first_diff(a, b):
    """First differing lines of two files, for the report."""
    out = subprocess.run(
        ["bash", "-c", f"diff '{a}' '{b}' | head -6"], capture_output=True, text=True
    )
    return out.stdout.rstrip()


def compare_exact(a, b):
    """Byte comparison with an order-only fallback."""
    if filecmp.cmp(a, b, shallow=False):
        return True, "identical"
    if _sorted_equal(a, b):
        return False, "same rows, different order"
    return False, "differ\n" + _first_diff(a, b)


def _normalize_sets(path):
    """Read a headerless TSV as strings, sorting items inside ';'-joined cells."""
    df = pl.read_csv(
        path, separator="\t", has_header=False, infer_schema_length=0, quote_char=None
    )
    return df.with_columns(
        pl.col(c).str.split(";").list.sort().list.join(";")
        for c in df.columns
        if df[c].str.contains(";").any()
    )


def compare_sets(a, b):
    """Row-for-row equality after normalizing ';'-joined columns."""
    x, y = _normalize_sets(a), _normalize_sets(b)
    if x.shape != y.shape:
        return False, f"shape {x.shape} vs {y.shape}"
    if x.equals(y):
        return True, "identical after set normalization"
    if x.sort(x.columns).equals(y.sort(y.columns)):
        return False, "same rows after set normalization, different order"
    bad = (x != y).select(pl.any_horizontal(pl.all()).sum()).item()
    return False, f"{bad} rows differ after set normalization"


def _read_z(path, has_header):
    """anchor + zscore as (str, f64), whitespace stripped; parquet read as-is."""
    if path.suffix == ".parquet":
        df = pl.read_parquet(path)
        return df.select(anchor=df.columns[0], z=pl.col(df.columns[1]).cast(pl.Float64))
    df = pl.read_csv(
        path, separator="\t", has_header=has_header, infer_schema_length=0,
        columns=[0, 1], quote_char=None,
    )
    a, z = df.columns
    return df.select(
        anchor=pl.col(a).str.strip_chars(),
        z=pl.col(z).str.strip_chars().cast(pl.Float64, strict=False),
    )


def compare_z(a, b, has_header, tol):
    """Same anchor set; z-scores within tol (nulls must match)."""
    x, y = _read_z(a, has_header), _read_z(b, has_header)
    if x.height != y.height:
        return False, f"rows {x.height} vs {y.height}"
    j = x.join(y, on="anchor", how="full", suffix="_new", coalesce=True)
    unmatched = j.filter(pl.col("z").is_null() & pl.col("z_new").is_null()).height
    missing = j.height - x.height
    if missing:
        return False, f"{missing} anchors not shared"
    null_mismatch = j.filter(pl.col("z").is_null() != pl.col("z_new").is_null()).height
    worst = j.select((pl.col("z") - pl.col("z_new")).abs().max()).item() or 0.0
    order = "same order" if x["anchor"].equals(y["anchor"]) else "different order"
    ok = null_mismatch == 0 and worst <= tol
    return ok, f"max |dz| {worst:.3g}, null mismatches {null_mismatch}, both-null {unmatched}, {order}"


def compare(a, b, mode, tol):
    """Dispatch on mode; returns (ok, message)."""
    if mode == "exact":
        return compare_exact(a, b)
    if mode == "rows":
        ok, msg = compare_exact(a, b)
        return ok or msg == "same rows, different order", msg
    if mode == "sets":
        return compare_sets(a, b)
    if mode == "zscore":
        return compare_z(a, b, has_header=False, tol=tol)
    if mode == "maxz":
        return compare_z(a, b, has_header=True, tol=tol)
    raise ValueError(f"Unknown mode: {mode}")


def expand(golden, new, g_rel, n_rel, key_regex=None):
    """Yield (name, golden_file, new_file) pairs; `{}` expands over golden matches."""
    if "{}" not in g_rel:
        yield g_rel, golden / g_rel, new / n_rel
        return
    head, tail = g_rel.split("{}")
    for g in sorted(golden.glob(g_rel.replace("{}", "*"))):
        rel = str(g.relative_to(golden))
        key = rel[len(head):len(rel) - len(tail)]
        new_key = re.match(key_regex, key).group(1) if key_regex else key
        yield g_rel.replace("{}", key), g, new / n_rel.replace("{}", new_key)


def main(argv=None):
    """Run every manifest check and exit non-zero if any fails."""
    args = build_parser().parse_args(argv)
    lines = [
        l.split("\t") for l in args.manifest.read_text().splitlines()
        if l.strip() and not l.startswith("#")
    ]
    failed = 0
    for g_rel, n_rel, mode, *key_regex in lines:
        g_rel = g_rel.replace("{prefix}", args.prefix)
        if args.only and not re.search(args.only, g_rel):
            continue
        n_ok = n_missing = n_bad = 0
        for name, g, n in expand(args.golden, args.new, g_rel, n_rel, *key_regex):
            if not n.exists():
                n_missing += 1
                continue
            ok, msg = compare(g, n, mode, args.tol)
            n_ok += ok
            n_bad += not ok
            if not ok or "{}" not in g_rel:
                print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}")
        if "{}" in g_rel:
            status = "FAIL" if n_bad else "PASS"
            print(f"[{status}] {g_rel}: {n_ok} pass, {n_bad} fail, {n_missing} not built")
        elif n_missing:
            print(f"[SKIP] {g_rel}: not built")
        failed += n_bad
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
