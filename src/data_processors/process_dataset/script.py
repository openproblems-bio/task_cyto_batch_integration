import sys
import anndata as ad
import openproblems as op

## VIASH START
par = {
    'input': 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/common_dataset.h5ad',
    'output_censored_left': 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_left.h5ad',
    'output_censored_right': 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_right.h5ad',
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

adata_unintegrated = adata.copy()

output_unintegrated = subset_h5ad_by_format(
    adata_unintegrated,
    config,
    "output_unintegrated"
)
print(f"output_unintegrated: {output_unintegrated}")

print(">> Creating split 1 data", flush=True)

output_censored_left = adata[(adata.obs.is_control>0) | (adata.obs.split==1)]

print("Grouping comparison:", flush=True)
print(output_censored_left.obs.groupby(["is_control", "split"]).size().to_dict())

output_censored_left = subset_h5ad_by_format(
    output_censored_left,
    config,
    "output_censored_left"
)
print(f"output_censored_left: {output_censored_left}")

print(">> Creating split 2 data", flush=True)

output_censored_right = adata[(adata.obs.is_control>0) | (adata.obs.split==2)]

print("Grouping comparison:", flush=True)
print(output_censored_right.obs.groupby(["is_control", "split"]).size().to_dict())

output_censored_right = subset_h5ad_by_format(
    output_censored_right,
    config,
    "output_censored_right"
)
print(f"output_censored_right: {output_censored_right}")

print(">> Writing data", flush=True)
output_unintegrated.write_h5ad(par["output_unintegrated"], compression="gzip")
output_censored_left.write_h5ad(par["output_censored_left"], compression="gzip")
output_censored_right.write_h5ad(par["output_censored_right"], compression="gzip")
