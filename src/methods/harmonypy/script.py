import anndata as ad
import harmonypy
import numpy as np

## VIASH START
par = {
    "input": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_split1.h5ad",
    "output": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/output_harmony_split1.h5ad",
}
meta = {"name": "harmonypy"}
## VIASH END

print("Reading and preparing input files", flush=True)
adata = ad.read_h5ad(par["input"])

# harmony can't handle integer batch labels
adata.obs["batch_str"] = adata.obs["batch"].astype(str)

markers_to_correct = adata.var[adata.var["to_correct"]].index.to_numpy()
markers_not_correct = adata.var[~adata.var["to_correct"]].index.to_numpy()

adata_to_correct = adata[:, markers_to_correct].copy()

print("Run harmony", flush=True)

# TODO numerical instability in kmeans causing problem with harmony.
# so adding a very small value to all entries to make sure there are no zeros
epsilon = 1e-20

out = harmonypy.run_harmony(
    data_mat=adata_to_correct.layers["preprocessed"] + epsilon,
    meta_data=adata_to_correct.obs,
    vars_use="batch_str",
)

# have to add in the uncorrected markers as well
uncorrected_data = adata[:, markers_not_correct].layers["preprocessed"]
out_matrix = np.concatenate([out.Z_corr, uncorrected_data], axis=1)
out_var_idx = np.concatenate([markers_to_correct, markers_not_correct])

# create new anndata
out_adata = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var.loc[out_var_idx][[]],
    layers={"integrated": out_matrix},
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "method_id": meta["name"],
        "parameters": {},
    },
)

# reorder var to match input
out_adata = out_adata[:, adata.var_names]

print("Write output AnnData to file", flush=True)

out_adata.write_h5ad(par["output"], compression="gzip")
