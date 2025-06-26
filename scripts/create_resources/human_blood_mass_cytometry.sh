#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DIR=resources_raw/VanGassen_dataset/
DATASET_ID=human_blood_mass_cytometry
OUTPUT_DIR=resources/datasets_raw/$DATASET_ID/

mkdir -p $OUTPUT_DIR

# create raw dataset files
python << HERE
import anndata as ad

adata = ad.read_h5ad("$RAW_DIR/VanGassen_dataset.h5ad")
adata.uns["dataset_id"] = "$DATASET_ID"
adata.uns["dataset_name"] = "Human Blood Mass Cytometry"

# override summary and description
adata.uns['dataset_summary'] = 'Mass cytometry dataset of whole blood from 2 human donors. ' \
  'For each donor, 2 samples are included: one unstimulated and one stimulated. ' \
  'For each of these, aliquotes of the same original sample were divided into 2 batches to allow the creation of sample-paired replicates for benchmarking purposes.'

adata.uns['dataset_description'] = 'CyTOF data (whole blood) from 2 human donors, each with unstimulated and stimulated samples. ' \
  'The stimulated samples were treated with Interferonα (IFNα) and Lipopolysaccharide (LPS). ' \
  'The mass cytometry cytometry panel includes 21 surface markers for cell phenotyping and 14 markers for functional analysis of the signaling responses. ' \
  'The data has been arcsinh transformed with cofactor 5 and pregated to only include events which fall either in the CD66-CD45+ or CD66+ gates.'

# rename values
adata.uns["parameter_som_xdim"] = int(adata.uns.pop("parameter_flowsom_xdim"))
adata.uns["parameter_som_ydim"] = int(adata.uns.pop("parameter_flowsom_ydim"))
adata.uns["parameter_num_clusters"] = int(adata.uns.pop("parameter_flowsom_nclus"))

adata.write_h5ad("$OUTPUT_DIR/common_dataset.h5ad", compression='gzip')
HERE

cat > $OUTPUT_DIR/state.yaml << HERE
id: $DATASET_ID
output_dataset: !file common_dataset.h5ad
HERE

# only run this if you have access to the openproblems-data bucket
aws s3 sync --profile op \
  resources/datasets_raw/$DATASET_ID \
  s3://openproblems-data/resources/task_cyto_batch_integration/datasets_raw/$DATASET_ID/ \
  --delete --dryrun

# # run the dataset processor
# cat > /tmp/params.yaml << 'HERE'
# input_states: s3://openproblems-data/resources/task_cyto_batch_integration/datasets_raw/**/state.yaml
# rename_keys: 'input:output_dataset'
# output_state: '$id/state.yaml'
# settings: '{"output_unintegrated": "$id/unintegrated.h5ad", "output_unintegrated_censored": "$id/unintegrated_censored.h5ad", "output_validation": "$id/validation.h5ad"}'
# publish_dir: s3://openproblems-data/resources/task_cyto_batch_integration/datasets/
# HERE

# tw launch https://github.com/openproblems-bio/task_cyto_batch_integration.git \
#   --revision build/main \
#   --pull-latest \
#   --main-script target/nextflow/workflows/process_datasets/main.nf \
#   --workspace 53907369739130 \
#   --params-file /tmp/params.yaml \
#   --entry-name auto \
#   --config common/nextflow_helpers/labels_tw.config \
#   --labels task_cyto_batch_integration,process_datasets
