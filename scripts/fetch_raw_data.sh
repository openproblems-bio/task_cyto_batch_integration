#!/bin/bash

set -e

echo "Syncing raw datasets from S3..."
aws s3 sync \
  s3://openproblems-data/resources/task_cyto_batch_integration/resources_raw/ \
  resources_raw
