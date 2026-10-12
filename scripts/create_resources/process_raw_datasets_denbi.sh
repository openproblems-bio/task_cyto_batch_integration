#!/bin/bash

# de.NBI variant of process_raw_datasets.sh: runs the process_datasets workflow on the
# BiBiGrid SLURM cluster and publishes straight to S3. The S3 credentials are the `tower`
# service user's, provided on the cluster via /vol/scratch/cluster-env.sh
# (AWS_SHARED_CREDENTIALS_FILE) — no credentials need to be passed at launch.

cat > /tmp/params.yaml << 'HERE'
input_states: s3://openproblems-data/resources/task_cyto_batch_integration/datasets_raw/**/state.yaml
rename_keys: 'input:output_dataset'
output_state: '$id/state.yaml'
settings: '{"output_unintegrated": "$id/unintegrated.h5ad", "output_censored_split1": "$id/censored_split1.h5ad", "output_censored_split2": "$id/censored_split2.h5ad"}'
publish_dir: s3://openproblems-data/resources/task_cyto_batch_integration/datasets/
HERE

tw launch https://github.com/openproblems-bio/task_cyto_batch_integration.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 3qstFmP9lNwdzutSNuJq7c \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config common/nextflow_helpers/labels_denbi.config \
  --labels task_cyto_batch_integration,process_datasets,denbi
