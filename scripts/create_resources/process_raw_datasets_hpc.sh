#!/bin/bash

# script to launch the process raw dataset workflow on slurm via seqera tower.

cat > /tmp/params.yaml << 'HERE'
input_states: s3://openproblems-data/resources/task_cyto_batch_integration/datasets_raw/**/state.yaml
rename_keys: 'input:output_dataset'
output_state: '$id/state.yaml'
settings: '{"output_unintegrated": "$id/unintegrated.h5ad", "output_censored_split1": "$id/censored_split1.h5ad", "output_censored_split2": "$id/censored_split2.h5ad"}'
publish_dir: /vast/scratch/users/putri.g/cytobenchmark/benchmark_out_hpc/datasets/
HERE

tw launch https://github.com/openproblems-bio/task_cyto_batch_integration.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/process_datasets/main.nf \
  --workspace 80689470953249 \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config scripts/labels_tw_wehi.config \
  --labels task_cyto_batch_integration,process_datasets
