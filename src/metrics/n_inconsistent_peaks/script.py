import anndata as ad
import scipy

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
  'input_validation': 'resources_test/.../validation.h5ad',
  'input_unintegrated': 'resources_test/.../unintegrated.h5ad',
  'input_integrated': 'resources_test/.../integrated.h5ad',
  'output': 'output.h5ad'
}
meta = {
  'name': 'n_inconsistent_peaks'
}
## VIASH END

print('Reading input files', flush=True)
input_validation = ad.read_h5ad(par['input_validation'])
input_unintegrated = ad.read_h5ad(par['input_unintegrated'])
input_integrated = ad.read_h5ad(par['input_integrated'])

print('Compute metrics', flush=True)
# metric_ids and metric_values can have length > 1
# but should be of equal length
uns_metric_ids = [ 'n_inconsistent_peaks' ]
uns_metric_values = [ 0.5 ]

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
    uns={
    'dataset_id': input_integrated.uns['dataset_id'],
    'method_id': input_integrated.uns['method_id'],
    "sample_ids": [],
    'metric_ids': uns_metric_ids,
    'metric_values': uns_metric_values
  }
)
output.write_h5ad(par['output'], compression='gzip')
