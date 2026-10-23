"""Generate the synthetic demo dataset, audit it, and refresh the README figures and sample report.

    python examples/run_demo.py
"""

import shutil
from pathlib import Path

from wellqc.pipeline import run_audit
from wellqc.synthetic import make_demo_dataset

ROOT = Path(__file__).resolve().parents[1]

if __name__ == "__main__":
    projects = make_demo_dataset(ROOT / "demo_data", seed=7)
    audit = run_audit(projects, ROOT / "demo_results")
    docs = ROOT / "docs"
    (docs / "figures").mkdir(parents=True, exist_ok=True)
    for f in audit.figures:
        shutil.copy(f, docs / "figures" / Path(f).name)
    report = audit.report_path.read_text()
    (docs / "sample_report.md").write_text(report)
    for name in ("well_master.csv", "issues.csv", "lineage.csv"):
        shutil.copy(audit.report_path.parent / name, docs / f"sample_{name}")
    print(f"Figures and sample report written to {docs}")
