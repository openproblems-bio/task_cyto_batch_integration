import sys
import anndata as ad
import openproblems as op

## VIASH START
par = {
    'input': 'resources_test/task_cyto_batch_integration/starter_file/common_dataset.h5ad',
    'validation_sample_names': [],
    'output_unintegrated': 'unintegrated.h5ad',
    'output_unintegrated_censored': 'unintegrated_censored.h5ad',
    'output_validation': 'validation.h5ad'
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

validation_names = par["validation_sample_names"] or []
is_validation = adata.obs["sample"].isin(validation_names)


print(">> Creating train data", flush=True)
output_unintegrated = subset_h5ad_by_format(
    adata[[not x for x in is_validation]],
    config,
    "output_unintegrated"
)
print(f"output_unintegrated: {output_unintegrated}")

print(">> Creating test data", flush=True)
output_unintegrated_censored = subset_h5ad_by_format(
    adata[[not x for x in is_validation]],
    config,
    "output_unintegrated_censored"
)
print(f"output_unintegrated_censored: {output_unintegrated_censored}")

print(">> Creating solution data", flush=True)
output_validation = subset_h5ad_by_format(
    adata[is_validation],
    config,
    "output_validation"
)
print(f"output_validation: {output_validation}")

print(">> Writing data", flush=True)
output_unintegrated.write_h5ad(par["output_unintegrated"], compression="gzip")
output_unintegrated_censored.write_h5ad(par["output_unintegrated_censored"], compression="gzip")
output_validation.write_h5ad(par["output_validation"], compression="gzip")
