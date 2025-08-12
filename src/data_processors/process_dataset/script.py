import sys
import anndata as ad
import openproblems as op

## VIASH START
par = {
    'input': 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/common_dataset.h5ad',
    'output_censored_split1': 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_split1.h5ad',
    'output_censored_split2': 'resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_split2.h5ad',
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

output_censored_split1 = adata[(adata.obs.is_control>0) | (adata.obs.split==1)]

print("Grouping comparison:", flush=True)
print(output_censored_split1.obs.groupby(["is_control", "split"]).size().to_dict())

output_censored_split1 = subset_h5ad_by_format(
    output_censored_split1,
    config,
    "output_censored_split1"
)
print(f"output_censored_split1: {output_censored_split1}")

print(">> Creating split 2 data", flush=True)

output_censored_split2 = adata[(adata.obs.is_control>0) | (adata.obs.split==2)]

print("Grouping comparison:", flush=True)
print(output_censored_split2.obs.groupby(["is_control", "split"]).size().to_dict())

output_censored_split2 = subset_h5ad_by_format(
    output_censored_split2,
    config,
    "output_censored_split2"
)
print(f"output_censored_split2: {output_censored_split2}")

print(">> Writing data", flush=True)
output_unintegrated.write_h5ad(par["output_unintegrated"], compression="gzip")
output_censored_split1.write_h5ad(par["output_censored_split1"], compression="gzip")
output_censored_split2.write_h5ad(par["output_censored_split2"], compression="gzip")
