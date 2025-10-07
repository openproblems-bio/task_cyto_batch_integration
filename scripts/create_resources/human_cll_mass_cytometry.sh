#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DIR=resources_raw/human_cll_mass_cytometry/
DATASET_ID=human_cll_mass_cytometry
OUTPUT_DIR=resources/datasets_raw/$DATASET_ID/

mkdir -p $OUTPUT_DIR

# create raw dataset files
python << HERE
import anndata as ad

adata = ad.read_h5ad("$RAW_DIR/human_cll_mass_cytometry.h5ad")

# make sure the output is compressed
adata.write_h5ad("$OUTPUT_DIR/common_dataset.h5ad", compression='gzip')
HERE

cat > $OUTPUT_DIR/state.yaml << HERE
id: $DATASET_ID
output_dataset: !file common_dataset.h5ad
HERE

# only run this if you have access to the openproblems-data bucket
aws s3 sync --profile op \
  resources/datasets_raw/human_cll_mass_cytometry \
  s3://openproblems-data/resources/task_cyto_batch_integration/datasets_raw/human_cll_mass_cytometry/ \
  --delete --dryrun
