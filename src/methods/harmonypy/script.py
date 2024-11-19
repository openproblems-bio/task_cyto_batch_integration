import anndata as ad
import harmonypy

## VIASH START
par = {
  'input': 'resources_test/task_cyto_batch_integration/starter_file/unintegrated_censored.h5ad',
  'output': 'output.h5ad'
}
meta = {
  'name': 'harmonypy'
}
## VIASH END

print('Reading input files', flush=True)
input = ad.read_h5ad(par['input'])

print('Run harmony', flush=True)
# harmony can't handle integer batch labels
input.obs["batch_str"] = input.obs["batch"].astype(str)
out = harmonypy.run_harmony(
  input.layers["preprocessed"],
  input.obs,
  "batch_str"
)

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  obs=input.obs[[]],
  var=input.var[[]],
  layers={
    "integrated": out.Z_corr.transpose()
  },
  uns={
    "dataset_id": input.uns["dataset_id"],
    "method_id": meta["name"],
    "parameters": {}
  }
)

output.write_h5ad(par['output'], compression='gzip')
