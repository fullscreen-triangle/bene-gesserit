"""
Master validation runner for the BPS membrane-computing series.

Usage (from phosphate/src/):
    python -m validation.run_all

Runs all six paper validation modules, collects results, and writes:
    phosphate/results/validation_results.json

Exit code 0 if all claims pass; 1 if any fail.
"""

import json
import sys
import os
import datetime
import platform
import numpy as np


class _NumpyEncoder(json.JSONEncoder):
    """Serialise numpy scalars and arrays to native Python types."""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.bool_):
            return bool(obj)
        return super().default(obj)

# Make sure the package is importable when run as a script from phosphate/src/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from validation import (                         # noqa: E402
    bilayer_substrate,
    state_encoder,
    aperture_array,
    energy_transducer,
    cascade_fabric,
    readout_interface,
)

PAPERS = [
    ("bilayer_substrate",  bilayer_substrate),
    ("state_encoder",      state_encoder),
    ("aperture_array",     aperture_array),
    ("energy_transducer",  energy_transducer),
    ("cascade_fabric",     cascade_fabric),
    ("readout_interface",  readout_interface),
]

RESULTS_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "results",
    "validation_results.json",
)


def run() -> dict:
    all_results = []
    per_paper   = {}

    for paper_key, module in PAPERS:
        print(f"\n{'-'*60}")
        print(f"  Paper: {paper_key}")
        print(f"{'-'*60}")

        results = module.run_all()
        per_paper[paper_key] = {
            "claims":  len(results),
            "passed":  sum(1 for r in results if r["passed"]),
            "failed":  sum(1 for r in results if not r["passed"]),
            "results": results,
        }
        all_results.extend(results)

        for r in results:
            status = "PASS" if r["passed"] else "FAIL"
            desc = r["description"][:62].encode("ascii", "replace").decode("ascii")
            print(f"  [{status}] {r['claim_id']:10s}  {desc}")

    total   = len(all_results)
    passed  = sum(1 for r in all_results if r["passed"])
    failed  = total - passed

    summary = {
        "timestamp":      datetime.datetime.utcnow().isoformat() + "Z",
        "python_version": platform.python_version(),
        "platform":       platform.platform(),
        "total_claims":   total,
        "passed":         passed,
        "failed":         failed,
        "pass_rate":      round(passed / total, 4) if total > 0 else 0.0,
    }

    output = {
        "summary":     summary,
        "papers":      per_paper,
        "all_results": all_results,
    }

    # ── Write JSON ────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
    with open(RESULTS_PATH, "w", encoding="utf-8") as fh:
        json.dump(output, fh, indent=2, ensure_ascii=False, cls=_NumpyEncoder)

    # ── Final summary ──────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"  TOTAL: {passed}/{total} passed  ({100*summary['pass_rate']:.1f}%)")
    print(f"  Results written to: {RESULTS_PATH}")
    print(f"{'='*60}\n")

    return output


def main() -> int:
    result = run()
    return 0 if result["summary"]["failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
