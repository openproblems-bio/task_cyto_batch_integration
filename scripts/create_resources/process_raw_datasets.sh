#!/bin/bash

cat > /tmp/params.yaml << 'HERE'
input_states: s3://openproblems-data/resources/task_cyto_batch_integration/datasets_raw/**/state.yaml
rename_keys: 'input:output_dataset'
output_state: '$id/state.yaml'
settings: '{"output_unintegrated": "$id/unintegrated.h5ad", "output_unintegrated_censored": "$id/unintegrated_censored.h5ad", "output_validation": "$id/validation.h5ad"}'
publish_dir: s3://openproblems-data/resources/task_cyto_batch_integration/datasets/
HERE

tw launch https://github.com/openproblems-bio/task_cyto_batch_integration.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config common/nextflow_helpers/labels_tw.config \
  --labels task_cyto_batch_integration,process_datasets
