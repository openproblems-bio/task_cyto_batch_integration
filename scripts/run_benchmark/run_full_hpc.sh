#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# generate a unique id
RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="/vast/scratch/users/putri.g/cytobenchmark/benchmark_out_hpc/results/${RUN_ID}"

# write the parameters to file
cat > /tmp/params.yaml << HERE
input_states: /vast/scratch/users/putri.g/cytobenchmark/benchmark_out_hpc/datasets/**/state.yaml
rename_keys: 'input_censored_split1:output_censored_split1;input_censored_split2:output_censored_split2;input_unintegrated:output_unintegrated'
output_state: "state.yaml"
settings: '{"metrics_include": ["ratio_inconsistent_peaks"], "methods_include": ["cytonorm_all_controls_to_goal"]}'
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_cyto_batch_integration.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 80689470953249 \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config scripts/labels_tw_wehi.config \
  --labels task_cyto_batch_integration,test_subset
