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

RAW_DATA=resources_test/task_cyto_batch_integration/starter_file

DATASET_DIR=resources_test/task_cyto_batch_integration/starter_file

mkdir -p $DATASET_DIR

# TODO: get original_dataset.h5ad from somewhere

# wget https://zenodo.org/records/13928969/files/ID1_Panel1_TP1.fcs?download=1 \
#   -O $DATASET_DIR/ID1_Panel1_TP1.fcs

# python << HERE
# import readfcs
# ad = readfcs.read("$DATASET_DIR/ID1_Panel1_TP1.fcs")
# ad.layers["transformed"] = ad.X
# del ad.X
# # todo: add other preprocessing steps to make sure the dataset is a common dataset
# ad.write_h5ad("$DATASET_DIR/common_dataset.h5ad")
# HERE

python << HERE
import anndata as ad

adata = ad.read_h5ad("resources_test/task_cyto_batch_integration/starter_file/original_dataset.h5ad")

channelsofinterest = ['UV379-A',
 'UV515-A',
 'UV610-A',
 'UV735-A',
 'V431-A',
 'V525-A',
 'V586-A',
 'V605-A',
 'V677-A',
 'V710-A',
 'V750-A',
 'V810-A',
 'B530-A',
 'B710-A',
 'YG586-A',
 'YG610-A',
 'YG670-A',
 'YG780-A',
 'R670-A',
 'R730-A']
adata.var.rename(columns={"n":"numeric_id"}, inplace=True)
marker_types = ["lineage" if chan in channelsofinterest else 'functional' for chan in adata.var["channel"]]
to_correct = [True if chan in channelsofinterest else False for chan in adata.var["channel"]]
adata.var["marker_type"] = marker_types
adata.var['to_correct'] = to_correct
adata.uns['dataset_id'] = 'XXXXX'
adata.uns['dataset_name'] = 'Summer School data'
adata.uns['dataset_summary'] = 'Draft data for cytometry batch integration benchmark'
adata.uns['dataset_description'] = '''
This is a draft dataset for the cytometry batch integration benchmark (Summer School). 
It contains only samples from one batch (Day1). 
Even though a preprocessed layer is available, it only contains arcsinh transformed data (not cleaned or compensated data).
'''
adata.uns['dataset_url'] = "https://saeyslab.sites.vib.be"
adata.uns['dataset_organism'] = "mus_musculus"
adata.uns['dataset_reference'] = "unpublished"

out_file = "resources_test/task_cyto_batch_integration/starter_file/common_dataset.h5ad"
adata.write_h5ad(out_file, compression="gzip")
HERE

# process dataset
viash run src/data_processors/process_dataset/config.vsh.yaml -- \
  --input $RAW_DATA/common_dataset.h5ad \
  --output_unintegrated $DATASET_DIR/unintegrated.h5ad \
  --output_unintegrated_censored $DATASET_DIR/unintegrated_censored.h5ad \
  --output_validation $DATASET_DIR/validation.h5ad

# # run one method
# viash run src/methods/logistic_regression/config.vsh.yaml -- \
#     --input_train $DATASET_DIR/cxg_mouse_pancreas_atlas/train.h5ad \
#     --input_test $DATASET_DIR/cxg_mouse_pancreas_atlas/test.h5ad \
#     --output $DATASET_DIR/cxg_mouse_pancreas_atlas/prediction.h5ad

# # run one metric
# viash run src/metrics/accuracy/config.vsh.yaml -- \
#     --input_prediction $DATASET_DIR/cxg_mouse_pancreas_atlas/prediction.h5ad \
#     --input_solution $DATASET_DIR/cxg_mouse_pancreas_atlas/solution.h5ad \
#     --output $DATASET_DIR/cxg_mouse_pancreas_atlas/score.h5ad

# # write manual state.yaml. this is not actually necessary but you never know it might be useful
# cat > $DATASET_DIR/cxg_mouse_pancreas_atlas/state.yaml << HERE
# id: cxg_mouse_pancreas_atlas
# train: !file train.h5ad
# test: !file test.h5ad
# solution: !file solution.h5ad
# prediction: !file prediction.h5ad
# score: !file score.h5ad
# HERE

# # only run this if you have access to the openproblems-data bucket
# aws s3 sync --profile op \
#   "$DATASET_DIR" s3://openproblems-data/resources_test/task_cyto_batch_integration \
#   --delete --dryrun
