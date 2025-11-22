import time

import anndata as ad
import numpy as np
import scvi
import torch
from scvi.external import cytovi

# from sklearn.cluster import KMeans
# from threadpoolctl import threadpool_limits

## VIASH START
par = {
    "input": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_split2.h5ad",
    "output": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/output_cytovi_split2.h5ad",
    "n_hidden": 128,
    "n_layers": 1,
}
meta = {"name": "cytovi"}
## VIASH END

# setting calculation to TF32 to speed up training
torch.backends.cuda.matmul.allow_tf32 = True

# increase num workers for data loading
scvi.settings.num_workers = 95

print("Reading and preparing input files", flush=True)
adata = ad.read_h5ad(par["input"])

adata.obs["batch_str"] = adata.obs["batch"].astype(str)
adata.obs["sample_key_str"] = adata.obs["sample"].astype(str)

markers_to_correct = adata.var[adata.var["to_correct"]].index.to_numpy()
markers_not_correct = adata.var[~adata.var["to_correct"]].index.to_numpy()

adata_to_correct = adata[:, markers_to_correct].copy()

print(
    f"Train CytoVI on {adata_to_correct.shape[0]} cells",
    flush=True,
)

cytovi.CYTOVI.setup_anndata(
    adata_to_correct,
    layer="preprocessed",
    batch_key="batch_str",
    sample_key="sample_key_str",
)

model = cytovi.CYTOVI(
    adata_to_correct, n_hidden=par["n_hidden"], n_layers=par["n_layers"]
)

print("Start training CytoVI model", flush=True)

start = time.time()
model.train(
    batch_size=8192,
)
end = time.time()
print(f"Training took {end - start:.2f} seconds", flush=True)

# get batch corrected data
print("Correcting data", flush=True)
corrected_data = model.get_normalized_expression(adata=adata_to_correct)

# have to add in the uncorrected markers as well
uncorrected_data = adata[:, markers_not_correct].layers["preprocessed"]

out_matrix = np.concatenate([corrected_data, uncorrected_data], axis=1)
out_var_idx = np.concatenate([corrected_data.columns, markers_not_correct])

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

# leave this here for debugging purposes
# run umap for quick check
# import scanpy as sc

# test_adata = ad.AnnData(
#     X=out_adata.layers["integrated"].toarray(),
#     obs=adata.obs,
#     var=adata.var,
# )
# test_adata = test_adata[:, markers_to_correct]
# sc.pp.neighbors(test_adata, use_rep="X")
# sc.tl.umap(test_adata)
# sc.pl.umap(test_adata, color="batch")

print("Write output AnnData to file", flush=True)

out_adata.write_h5ad(par["output"], compression="gzip")
