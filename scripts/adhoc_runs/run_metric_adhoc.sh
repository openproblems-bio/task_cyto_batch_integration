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
    --output-dir /Users/putri.g/Documents/cytobenchmark/analysis/data/abundance_preservation/2026-07-30 \
    --skip-existing \
    --datasets mouse_spleen_flow_cytometry lille_spectral_flow_cytometry \
    --methods cycombine_all_controls_to_goal cycombine_all_controls_to_mid \
        cycombine_no_controls_to_goal cycombine_no_controls_to_mid \
        cycombine_one_control_to_goal cycombine_one_control_to_mid \
        cytonorm_all_controls_to_mid cytonorm_no_controls_to_mid cytonorm_one_control_to_mid
