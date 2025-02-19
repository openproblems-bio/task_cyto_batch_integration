#!/bin/bash

set -e

common/scripts/create_component \
  --name no_integration \
  --language python \
  --type control_method
