#!/usr/bin/env python3
"""
Run a Python metric ad-hoc for all (dataset, method) combinations found in a
pipeline output directory, without re-running the full pipeline.

The script discovers method outputs under:
    <input-dir>/<dataset>/method_out/<method>_split1.h5ad
    <input-dir>/<dataset>/method_out/<method>_split2.h5ad

and writes metric outputs to:
    <output-dir>/<dataset>/<metric>_<method>.h5ad

It works by replacing the ## VIASH START...END block in the metric script with
the correct par/meta values and running the result in a subprocess. This means
any Python metric in src/metrics/ can be run ad-hoc with this script.

Usage:
    python scripts/run_metric_adhoc.py \\
        --metric functional_marker_preservation \\
        --input-dir /path/to/run_2026-04-26_21-34-33 \\
        --datasets-dir /path/to/datasets \\
        --output-dir /path/to/analysis_2026-05-13/functional_marker_preservation

    # Restrict to specific datasets or methods:
    python scripts/run_metric_adhoc.py \\
        --metric functional_marker_preservation \\
        --input-dir /path/to/run_2026-04-26_21-34-33 \\
        --datasets-dir /path/to/datasets \\
        --output-dir /path/to/output \\
        --datasets mouse_spleen_flow_cytometry \\
        --methods harmonypy combat

    # Skip already-computed outputs:
    python scripts/run_metric_adhoc.py ... --skip-existing
"""

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.resolve()


def find_method_splits(input_dir: Path, dataset: str, method: str):
    """
    Locate the split1 and split2 output files for a (dataset, method) pair.

    Returns (split1_path, split2_path) or None if either file is missing.
    """
    method_out_dir = input_dir / dataset / "method_out"
    split1 = method_out_dir / f"{method}_split1.h5ad"
    split2 = method_out_dir / f"{method}_split2.h5ad"

    if not split1.exists():
        print(f"    WARNING: split1 not found: {split1}", flush=True)
        return None
    if not split2.exists():
        print(f"    WARNING: split2 not found: {split2}", flush=True)
        return None

    return split1, split2


def run_metric(metric_name: str, par: dict) -> None:
    """
    Inject par/meta into the metric script and execute it in a subprocess.

    The ## VIASH START...END block in the script is replaced with the provided
    par and meta values. The metric's own directory is added to sys.path inside
    the injected block so that `import helper` resolves correctly regardless of
    where this runner is called from.
    """
    script_path = REPO_ROOT / "src" / "metrics" / metric_name / "script.py"
    if not script_path.exists():
        raise FileNotFoundError(
            f"Metric script not found: {script_path}\n"
            f"Check that '{metric_name}' is a valid metric under src/metrics/."
        )

    metric_dir = script_path.parent
    meta = {
        "name": metric_name,
        "resources_dir": str(REPO_ROOT / "src" / "utils"),
    }

    # The injected block sets par/meta and ensures the metric's helper.py is
    # importable. sys is already imported at the top of every metric script, so
    # we can call sys.path.insert directly without re-importing it.
    injected_block = (
        f"## VIASH START\n"
        f"par = {par!r}\n"
        f"meta = {meta!r}\n"
        f"sys.path.insert(0, {str(metric_dir)!r})\n"
        f"## VIASH END"
    )

    original = script_path.read_text()
    patched = re.sub(
        r"## VIASH START.*?## VIASH END",
        injected_block,
        original,
        flags=re.DOTALL,
    )

    # Write the patched script to a temp file and run it as a subprocess.
    # Using a temp file (rather than exec/runpy) gives clean isolation: each
    # metric run gets its own interpreter state and import namespace.
    with tempfile.NamedTemporaryFile(
        suffix=".py", mode="w", delete=False, dir=metric_dir
    ) as f:
        f.write(patched)
        temp_path = Path(f.name)

    try:
        subprocess.run([sys.executable, str(temp_path)], check=True)
    finally:
        temp_path.unlink(missing_ok=True)


def discover_methods(method_out_dir: Path) -> list[str]:
    """Return all method names found as <method>_split1.h5ad in method_out_dir."""
    return sorted(
        p.stem.replace("_split1", "") for p in method_out_dir.glob("*_split1.h5ad")
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Run a Python metric ad-hoc for all (dataset, method) combinations "
            "in a pipeline output directory."
        )
    )
    parser.add_argument(
        "--metric",
        required=True,
        help="Metric name matching a directory under src/metrics/ (e.g. functional_marker_preservation).",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Pipeline output directory containing <dataset>/method_out/ subdirectories.",
    )
    parser.add_argument(
        "--datasets-dir",
        type=Path,
        required=True,
        help=(
            "Base directory containing unintegrated data, "
            "structured as <datasets-dir>/<dataset>/unintegrated.h5ad."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory to write metric output files into, as <output-dir>/<dataset>/<metric>_<method>.h5ad.",
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=None,
        metavar="DATASET",
        help="Restrict to these datasets (default: all found in input-dir).",
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        default=None,
        metavar="METHOD",
        help="Restrict to these methods (default: all found for each dataset).",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip (dataset, method) combinations where the output file already exists.",
    )
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    datasets_dir = args.datasets_dir.resolve()
    output_dir = args.output_dir.resolve()

    datasets = args.datasets or sorted(
        d.name for d in input_dir.iterdir() if d.is_dir()
    )

    for dataset in datasets:
        print(f"\nDataset: {dataset}", flush=True)

        method_out_dir = input_dir / dataset / "method_out"
        if not method_out_dir.exists():
            print("  No method_out directory found, skipping.", flush=True)
            continue

        unintegrated = datasets_dir / dataset / "unintegrated.h5ad"
        if not unintegrated.exists():
            print(
                f"  Unintegrated file not found: {unintegrated}, skipping.", flush=True
            )
            continue

        methods = args.methods or discover_methods(method_out_dir)
        if not methods:
            print(
                f"  No method outputs found in {method_out_dir}, skipping.", flush=True
            )
            continue

        for method in methods:
            print(f"  Method: {method}", flush=True)

            output_path = output_dir / dataset / f"{args.metric}_{method}.h5ad"
            if args.skip_existing and output_path.exists():
                print(f"    Output exists, skipping: {output_path}", flush=True)
                continue

            splits = find_method_splits(input_dir, dataset, method)
            if splits is None:
                continue

            split1, split2 = splits
            output_path.parent.mkdir(parents=True, exist_ok=True)

            par = {
                "input_unintegrated": str(unintegrated),
                "input_integrated_split1": str(split1),
                "input_integrated_split2": str(split2),
                "output": str(output_path),
            }

            print(f"    Output: {output_path}", flush=True)
            run_metric(args.metric, par)
            print("    Done.", flush=True)

    print("\nAll done.", flush=True)


if __name__ == "__main__":
    main()
