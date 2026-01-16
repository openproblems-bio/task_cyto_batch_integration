#!/bin/bash

# script to run each method as a job, which allows me to pull all the images to the HPC beforehand

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

METHODS=(
    "batchadjust_all_controls"
    "batchadjust_one_control"
    "combat"
    "cycombine_all_controls_to_goal"
    "cycombine_all_controls_to_mid"
    "cycombine_no_controls_to_goal"
    "cycombine_no_controls_to_mid"
    "cycombine_one_control_to_goal"
    "cycombine_one_control_to_mid"
    "cytonorm_all_controls_to_goal"
    "cytonorm_all_controls_to_mid"
    "cytonorm_no_controls_to_goal"
    "cytonorm_no_controls_to_mid"
    "cytonorm_one_control_to_goal"
    "cytonorm_one_control_to_mid"
    "cytovi"
    "gaussnorm"
    "harmonypy"
    "limma_remove_batch_effect"
    "rpca_to_goal"
    "rpca_to_mid"
    "no_integration"
    "perfect_integration"
    "shuffle_integration"
    "shuffle_integration_by_batch"
    "shuffle_integration_by_cell_type"
)

# Loop through each method
for method in "${METHODS[@]}"; do
  echo "============================================"
  echo "Submitting job for method: $method"
  echo "============================================"
  
  # generate a unique id
  RUN_ID="warmup_${method}_$(date +%Y-%m-%d_%H-%M-%S)"
  publish_dir="/vast/scratch/users/putri.g/cytobenchmark/benchmark_out_hpc/results/${RUN_ID}"

  # write the parameters to file
  cat > /tmp/params_${method}.yaml << HERE
input_states: /vast/scratch/users/putri.g/cytobenchmark/benchmark_out_hpc/datasets/**/state.yaml
rename_keys: 'input_censored_split1:output_censored_split1;input_censored_split2:output_censored_split2;input_unintegrated:output_unintegrated'
output_state: "state.yaml"
settings: '{"metrics_include": ["emd"], "methods_include": ["${method}"]}'
publish_dir: "$publish_dir"
HERE

  tw launch https://github.com/openproblems-bio/task_cyto_batch_integration.git \
    --revision build/setup_run_hpc \
    --pull-latest \
    --main-script target/nextflow/workflows/run_benchmark/main.nf \
    --workspace 80689470953249 \
    --params-file /tmp/params_${method}.yaml \
    --entry-name auto \
    --config scripts/labels_tw_wehi.config \
    --labels task_cyto_batch_integration,${method}
  
  echo "Job submitted for $method with RUN_ID: $RUN_ID"
  echo ""
  
  # Optional: Add a small delay to avoid overwhelming Tower
  sleep 5
done

echo "All method jobs submitted!"

# repeat for metrics
# do one for metrics as well later.
# emd and lisi has been used
# METRICS=(
#     "average_batch_r2"
#     "flowsom_mapping_similarity"
#     "ratio_inconsistent_peaks"
# )
# for metric in "${METRICS[@]}"; do
#   echo "============================================"
#   echo "Submitting job for metric: $metric"
#   echo "============================================"
  
#   # generate a unique id
#   RUN_ID="warmup_${metric}_$(date +%Y-%m-%d_%H-%M-%S)"
#   publish_dir="/vast/scratch/users/putri.g/cytobenchmark/benchmark_out_hpc/results/${RUN_ID}"

#   # write the parameters to file
#   cat > /tmp/params_${metric}.yaml << HERE
# input_states: /vast/scratch/users/putri.g/cytobenchmark/benchmark_out_hpc/datasets/**/state.yaml
# rename_keys: 'input_censored_split1:output_censored_split1;input_censored_split2:output_censored_split2;input_unintegrated:output_unintegrated'
# output_state: "state.yaml"
# settings: '{"metrics_include": ["${metric}"], "methods_include": ["cytonorm_all_controls_to_goal"]}'
# publish_dir: "$publish_dir"
# HERE

#   tw launch https://github.com/openproblems-bio/task_cyto_batch_integration.git \
#     --revision build/setup_run_hpc \
#     --pull-latest \
#     --main-script target/nextflow/workflows/run_benchmark/main.nf \
#     --workspace 80689470953249 \
#     --params-file /tmp/params_${metric}.yaml \
#     --entry-name auto \
#     --config scripts/labels_tw_wehi.config \
#     --labels task_cyto_batch_integration,${metric}
  
#   echo "Job submitted for $metric with RUN_ID: $RUN_ID"
#   echo ""
  
#   # Optional: Add a small delay to avoid overwhelming Tower
#   sleep 5
# done

# echo "All metrics jobs submitted!"