import sys
import anndata as ad
import openproblems as op

## VIASH START
par = {
    'input': 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/common_dataset.h5ad',
    'output_unintegrated': 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/unintegrated.h5ad',
    'output_unintegrated_censored': 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/unintegrated_censored.h5ad',
    'output_validation': 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/validation.h5ad'
}
meta = {
    'resources_dir': 'target/executable/data_processors/process_dataset',
    'config': 'target/executable/data_processors/process_dataset/.config.vsh.yaml'
}
## VIASH END

# import helper functions
sys.path.append(meta['resources_dir'])
from subset_h5ad_by_format import subset_h5ad_by_format

config = op.project.read_viash_config(meta["config"])


print(">> Load data", flush=True)
adata = ad.read_h5ad(par["input"])
print("input:", adata)

print(">> Creating unintegrated data", flush=True)

adata_unintegrated = adata[adata.obs.is_validation==False]

output_unintegrated = subset_h5ad_by_format(
    adata_unintegrated,
    config,
    "output_unintegrated"
)
print(f"output_unintegrated: {output_unintegrated}")

print(">> Creating test data", flush=True)
output_unintegrated_censored = subset_h5ad_by_format(
    adata_unintegrated,
    config,
    "output_unintegrated_censored"
)
print(f"output_unintegrated_censored: {output_unintegrated_censored}")

print(">> Creating validation data", flush=True)

adata_validation = adata[adata.obs.is_validation==True]

output_validation = subset_h5ad_by_format(
    adata_validation,
    config,
    "output_validation"
)
print(f"output_validation: {output_validation}")

print(">> Writing data", flush=True)
output_unintegrated.write_h5ad(par["output_unintegrated"], compression="gzip")
output_unintegrated_censored.write_h5ad(par["output_unintegrated_censored"], compression="gzip")
output_validation.write_h5ad(par["output_validation"], compression="gzip")
