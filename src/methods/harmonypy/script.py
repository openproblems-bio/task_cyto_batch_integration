import anndata as ad
import harmonypy

## VIASH START
par = {
    "input": "resources_test/task_cyto_batch_integration/starter_file/unintegrated_censored.h5ad",
    "output": "output.h5ad",
}
meta = {"name": "harmonypy"}
## VIASH END

print("Reading and preparing input files", flush=True)
adata = ad.read_h5ad(par["input"])

adata.obs["batch_str"] = adata.obs["batch"].astype(str)

markers_to_correct = adata.var[adata.var["to_correct"]].index.to_numpy()

adata_to_correct = adata[:, markers_to_correct]

print("Run harmony", flush=True)
# harmony can't handle integer batch labels

out = harmonypy.run_harmony(
    data_mat=adata_to_correct.layers["preprocessed"],
    meta_data=adata_to_correct.obs,
    vars_use="batch_str",
)

# create new anndata
out_adata = ad.AnnData(
    obs=adata.obs[[]],
    var=adata_to_correct.var[[]],
    layers={"integrated": out.Z_corr.transpose()},
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "method_id": meta["name"],
        "parameters": {},
    },
)

print("Write output AnnData to file", flush=True)

out_adata.write_h5ad(par["output"], compression="gzip")
