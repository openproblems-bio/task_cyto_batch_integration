#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# # remove this when you have implemented the script
# echo "TODO: replace the commands in this script with the sequence of components that you need to run to generate test_resources."
# echo "  Inside this script, you will need to place commands to generate example files for each of the 'src/api/file_*.yaml' files."
# exit 1

set -e

RAW_DATA=resources_test/task_cyto_batch_integration/leomazzi_cyto_spleen_subset

DATASET_DIR=resources_test/task_cyto_batch_integration/leomazzi_cyto_spleen_subset

mkdir -p $DATASET_DIR

python << HERE
import anndata as ad

adata = ad.read_h5ad("$RAW_DATA/common_dataset.h5ad")
adata.uns["dataset_id"] = "leomazzi_cyto_spleen_subset"
adata.uns["dataset_name"] = "Leomazzi Spleen Cytometry Subset"
adata.write_h5ad("$DATASET_DIR/common_dataset.h5ad")
HERE

# process dataset
viash run src/data_processors/process_dataset/config.vsh.yaml -- \
  --input $RAW_DATA/common_dataset.h5ad \
  --output_unintegrated $DATASET_DIR/unintegrated.h5ad \
  --output_unintegrated_censored $DATASET_DIR/unintegrated_censored.h5ad \
  --output_validation $DATASET_DIR/validation.h5ad

# run one method
viash run src/methods/harmonypy/config.vsh.yaml -- \
  --input $DATASET_DIR/unintegrated.h5ad \
  --output $DATASET_DIR/integrated.h5ad

# run one metric
viash run src/metrics/emd/config.vsh.yaml -- \
    --input_validation $DATASET_DIR/validation.h5ad \
    --input_unintegrated $DATASET_DIR/unintegrated.h5ad \
    --input_integrated $DATASET_DIR/integrated.h5ad \
    --output $DATASET_DIR/score.h5ad

# write manual state.yaml
cat > $DATASET_DIR/state.yaml << HERE
id: leomazzi_cyto_spleen_subset
integrated: !file integrated.h5ad
unintegrated: !file unintegrated.h5ad
unintegrated_censored: !file unintegrated_censored.h5ad
validation: !file validation.h5ad
score: !file score.h5ad
HERE

# only run this if you have access to the openproblems-data bucket
aws s3 sync --profile op \
  resources_test/task_cyto_batch_integration/leomazzi_cyto_spleen_subset \
  s3://openproblems-data/resources_test/task_cyto_batch_integration/leomazzi_cyto_spleen_subset \
  --delete --dryrun
