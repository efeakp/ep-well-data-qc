"""Command line interface.

    wellqc audit PROJECT_DIR [PROJECT_DIR ...] --out results/
    wellqc demo --out demo_results/
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _summary(a) -> str:
    sev = pd.Series([i.severity for i in a.issues]).value_counts()
    return (f"{len(a.inventory)} files, {len(a.master)} physical wells ({len(a.wells)} well records), "
            f"{sev.get('critical', 0)} critical / {sev.get('major', 0)} major / {sev.get('minor', 0)} minor issues\n"
            f"Report: {a.report_path}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="wellqc", description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    au = sub.add_parser("audit", help="audit one or more project folders (each folder = one legacy project)")
    au.add_argument("projects", nargs="+", type=Path)
    au.add_argument("--out", type=Path, default=Path("wellqc_results"))
    au.add_argument("--depth-unit", choices=["m", "ft"], default="m")
    au.add_argument("--uwi-prefix", default="W")
    au.add_argument("--no-predict", action="store_true")
    au.add_argument("--no-figures", action="store_true")
    de = sub.add_parser("demo", help="generate the synthetic demo dataset and audit it")
    de.add_argument("--out", type=Path, default=Path("demo_results"))
    de.add_argument("--data", type=Path, default=Path("demo_data"))
    de.add_argument("--seed", type=int, default=7)
    args = ap.parse_args(argv)

    from .pipeline import run_audit

    if args.cmd == "audit":
        projects = {p.name: p for p in args.projects}
        a = run_audit(projects, args.out, args.depth_unit, predict=not args.no_predict,
                      figures=not args.no_figures, uwi_prefix=args.uwi_prefix)
    else:
        from .synthetic import make_demo_dataset

        projects = make_demo_dataset(args.data, seed=args.seed)
        a = run_audit(projects, args.out)
    print(_summary(a))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
