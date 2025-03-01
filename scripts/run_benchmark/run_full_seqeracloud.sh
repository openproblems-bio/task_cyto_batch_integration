#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# generate a unique id
RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="s3://openproblems-data/resources/task_cyto_batch_integration/results/${RUN_ID}"

# write the parameters to file
cat > /tmp/params.yaml << HERE
# input_states: s3://openproblems-data/resources/task_cyto_batch_integration/datasets/**/state.yaml
input_states: s3://openproblems-data/resources_test/task_cyto_batch_integration/**/state.yaml
rename_keys: 'input_unintegrated_censored:unintegrated_censored;input_unintegrated:unintegrated;input_validation:validation'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_cyto_batch_integration.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 7gRyww9YNGb0c6BUBtLhDP \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config common/nextflow_helpers/labels_tw.config \
  --labels task_cyto_batch_integration,full
