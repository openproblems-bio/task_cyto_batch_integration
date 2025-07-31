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

RAW_DATA=resources_raw/Leomazzi_dataset
DATASET_ID=mouse_spleen_flow_cytometry_subset
DATASET_DIR=resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset

mkdir -p $DATASET_DIR

python << HERE
import anndata as ad

adata = ad.read_h5ad("$RAW_DATA/common_dataset.h5ad")

# determine split from is_validation and is_control
# if is_control >= 1 → 0
# else if not is_validation → 1
# else → 2
adata.obs["split"] = 1
adata.obs.loc[adata.obs["is_validation"], "split"] = 2
adata.obs.loc[adata.obs["is_control"] >= 1, "split"] = 0

# override dataset_id and dataset_name
adata.uns["dataset_id"] = "$DATASET_ID"
adata.uns["dataset_name"] = "Mouse Spleen Flow Cytometry Subset"

# rename values
adata.uns["parameter_som_xdim"] = int(adata.uns.pop("parameter_flowsom_xdim"))
adata.uns["parameter_som_ydim"] = int(adata.uns.pop("parameter_flowsom_ydim"))
adata.uns["parameter_num_clusters"] = int(adata.uns.pop("parameter_flowsom_nclus"))

adata.write_h5ad("$DATASET_DIR/common_dataset.h5ad", compression='gzip')
HERE

# process dataset
viash run src/data_processors/process_dataset/config.vsh.yaml -- \
  --input $DATASET_DIR/common_dataset.h5ad \
  --output_unintegrated $DATASET_DIR/unintegrated.h5ad \
  --output_censored_left $DATASET_DIR/censored_left.h5ad \
  --output_censored_right $DATASET_DIR/censored_right.h5ad

# run one method
viash run src/methods/harmonypy/config.vsh.yaml -- \
  --input $DATASET_DIR/censored_left.h5ad \
  --output $DATASET_DIR/integrated_left.h5ad

# run one method
viash run src/methods/harmonypy/config.vsh.yaml -- \
  --input $DATASET_DIR/censored_right.h5ad \
  --output $DATASET_DIR/integrated_right.h5ad

# run one metric
viash run src/metrics/emd/config.vsh.yaml -- \
    --input_unintegrated $DATASET_DIR/unintegrated.h5ad \
    --input_integrated_left $DATASET_DIR/integrated_left.h5ad \
    --input_integrated_right $DATASET_DIR/integrated_right.h5ad \
    --output $DATASET_DIR/score.h5ad

# write manual state.yaml
cat > $DATASET_DIR/state.yaml << HERE
id: $DATASET_ID
unintegrated: !file unintegrated.h5ad
censored_left: !file censored_left.h5ad
censored_right: !file censored_right.h5ad
integrated_left: !file integrated_left.h5ad
integrated_right: !file integrated_right.h5ad
score: !file score.h5ad
HERE

# only run this if you have access to the openproblems-data bucket
aws s3 sync --profile op \
  resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset \
  s3://openproblems-data/resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset \
  --delete --dryrun
