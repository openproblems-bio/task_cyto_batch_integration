import anndata as ad

## VIASH START
# Note: this section is auto-generated by viash at runtime. To edit it, make changes
# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
par = {
  "input_validation": "resources_test/task_cyto_batch_integration/starter_file/validation.h5ad",
  "output": "resources_test/task_cyto_batch_integration/starter_file/output.h5ad"
}
meta = {
  'name': 'perfect_integration'
}
## VIASH END

print("Reading input files", flush=True)
adata = ad.read_h5ad(par['input_validation'])

print("Creating integrated data", flush=True)

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    layers={"integrated": adata.layers['preprocessed']},
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "method_id": meta["name"],
        "parameters": {},
    },
)
output.write_h5ad(par['output'], compression='gzip')
