import anndata as ad
import mnnpy
import numpy as np

## VIASH START
par = {
    "input": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_split1.h5ad",
    "output": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/output_mnnpy.h5ad",
}
meta = {"name": "mnnpy"}
## VIASH END

print("Read input", flush=True)
adata = ad.read_h5ad(par["input"])

adata.X = adata.layers["preprocessed"]

# convert batch to category as otherwise mnnpy won't work..
adata.obs["batch_cat"] = adata.obs["batch"].astype(str)

markers_to_correct = adata.var[adata.var["to_correct"]].index.to_numpy()
markers_not_correct = adata.var[~adata.var["to_correct"]].index.to_numpy()

adata_to_correct = adata[:, markers_to_correct].copy()


print("Run mnn", flush=True)
split = []
batch_categories = adata_to_correct.obs["batch_cat"].unique().tolist()

for i in batch_categories:
    split.append(adata_to_correct[adata_to_correct.obs["batch_cat"] == i].copy())

corrected, _, _ = mnnpy.mnn_correct(
    *split, batch_key="batch", batch_categories=batch_categories, index_unique=None
)

# have to add in the uncorrected markers as well
uncorrected_data = adata[:, markers_not_correct].layers["preprocessed"]

out_matrix = np.concatenate([corrected.X, uncorrected_data], axis=1)
out_var_idx = np.concatenate([corrected.var.index, markers_not_correct])

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

# run umap for quick check
# import scanpy as sc
# test_adata = out_adata.copy()
# test_adata = test_adata[:, markers_to_correct]
# test_adata.X = test_adata.layers["integrated"]
# test_adata.obs = adata.obs
# sc.pp.neighbors(test_adata, use_rep="X")
# sc.tl.umap(test_adata)
# sc.pl.umap(test_adata, color="batch")


print("Store outputs", flush=True)
out_adata.write_h5ad(par["output"], compression="gzip")
