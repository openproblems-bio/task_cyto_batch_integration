#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/" && pwd)"

# python "$SCRIPT_DIR/run_metric_adhoc.py" \
#     --metric functional_marker_preservation \
#     --input-dir /Users/putri.g/Documents/cytobenchmark/analysis/runs/run_2026-04-26_21-34-33 \
#     --datasets-dir /Users/putri.g/Documents/cytobenchmark/datasets \
#     --output-dir /Users/putri.g/Documents/cytobenchmark/analysis/data/fmp/2026-07-15 \
#     --skip-existing
#     --datasets lille_spectral_flow_cytometry mouse_spleen_flow_cytometry

Rscript "$SCRIPT_DIR/run_metric_adhoc.R" \
    --metric abundance_preservation \
    --input-dir /Users/putri.g/Documents/cytobenchmark/analysis/runs/run_2026-04-26_21-34-33 \
    --datasets-dir /Users/putri.g/Documents/cytobenchmark/datasets \
    --output-dir /Users/putri.g/Documents/cytobenchmark/analysis/data/abundance_preservation/2026-08-25 \
    --skip-existing \
    --datasets mouse_spleen_flow_cytometry lille_spectral_flow_cytometry 