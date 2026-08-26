#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

resources_test_s3=s3://openproblems-data/resources_test/task_cyto_batch_integration

# generate a unique id
RUN_ID="testrun_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="/vol/scratch/results/task_cyto_batch_integration/${RUN_ID}"

# write the parameters to file
cat > /tmp/params.yaml << HERE
id: mouse_spleen_flow_cytometry_subset
input_unintegrated: $resources_test_s3/mouse_spleen_flow_cytometry_subset/unintegrated.h5ad
input_censored_split1: $resources_test_s3/mouse_spleen_flow_cytometry_subset/censored_split1.h5ad
input_censored_split2: $resources_test_s3/mouse_spleen_flow_cytometry_subset/censored_split2.h5ad
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/task_cyto_batch_integration.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 3qstFmP9lNwdzutSNuJq7c \
  --params-file /tmp/params.yaml \
  --config common/nextflow_helpers/labels_denbi.config \
  --labels task_cyto_batch_integration,test,denbi
