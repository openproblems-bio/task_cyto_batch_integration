import anndata as ad

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

print('Preprocess data', flush=True)
# ... preprocessing ...

print('Train model', flush=True)
# ... train model ...

print('Generate predictions', flush=True)
# ... generate predictions ...

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  
)
output.write_h5ad(par['output'], compression='gzip')
