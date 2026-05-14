#!/bin/bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

python "$REPO_ROOT/scripts/run_metric_adhoc.py" \
    --metric functional_marker_preservation \
    --input-dir /Users/putri.g/Documents/cytobenchmark/analysis/run_2026-04-26_21-34-33 \
    --datasets-dir /Users/putri.g/Documents/cytobenchmark/datasets \
    --output-dir /Users/putri.g/Documents/cytobenchmark/analysis/analysis_2026-05-13/functional_marker_preservation \
    --skip-existing
