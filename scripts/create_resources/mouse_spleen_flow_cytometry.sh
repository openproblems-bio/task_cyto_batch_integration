#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DIR=resources_raw/Leomazzi_dataset/
DATASET_ID=mouse_spleen_flow_cytometry
OUTPUT_DIR=resources/datasets_raw/$DATASET_ID/

mkdir -p $OUTPUT_DIR

# create raw dataset files
python << HERE
import anndata as ad

adata = ad.read_h5ad("$RAW_DIR/Leomazzi_dataset.h5ad")

# determine split from is_validation and is_control
# if is_control >= 1 → 0
# else if not is_validation → 1
# else → 2
adata.obs["split"] = 1
adata.obs.loc[adata.obs["is_validation"], "split"] = 2
adata.obs.loc[adata.obs["is_control"] >= 1, "split"] = 0

# add goal_batch
adata.uns["goal_batch"] = 1

# override dataset_id and dataset_name
adata.uns["dataset_id"] = "$DATASET_ID"
adata.uns["dataset_name"] = "Mouse Spleen Flow Cytometry"

# rename values
adata.uns["parameter_som_xdim"] = int(adata.uns.pop("parameter_flowsom_xdim"))
adata.uns["parameter_som_ydim"] = int(adata.uns.pop("parameter_flowsom_ydim"))
adata.uns["parameter_num_clusters"] = int(adata.uns.pop("parameter_flowsom_nclus"))

# make sure the output is compressed
adata.write_h5ad("$OUTPUT_DIR/common_dataset.h5ad", compression='gzip')
HERE

cat > $OUTPUT_DIR/state.yaml << HERE
id: $DATASET_ID
output_dataset: !file common_dataset.h5ad
HERE

# only run this if you have access to the openproblems-data bucket
aws s3 sync --profile op \
  resources/datasets_raw/mouse_spleen_flow_cytometry \
  s3://openproblems-data/resources/task_cyto_batch_integration/datasets_raw/mouse_spleen_flow_cytometry/ \
  --delete --dryrun
