import anndata as ad
import numpy as np
from scvi.external import cytovi

## VIASH START
par = {
    "input": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_split2.h5ad",
    "output": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/output_cytovi_split2.h5ad",
    "n_hidden": 128,
    "n_layers": 1,
}
meta = {"name": "cytovi"}
## VIASH END

print("Reading and preparing input files", flush=True)
adata = ad.read_h5ad(par["input"])

adata.obs["batch_str"] = adata.obs["batch"].astype(str)

markers_to_correct = adata.var[adata.var["to_correct"]].index.to_numpy()
markers_not_correct = adata.var[~adata.var["to_correct"]].index.to_numpy()

adata_to_correct = adata[:, markers_to_correct].copy()

# scale data
cytovi.scale(
    adata=adata_to_correct, transformed_layer_key="preprocessed", batch_key="batch_str"
)

print("Run CytoVI", flush=True)

cytovi.CYTOVI.setup_anndata(adata_to_correct, layer="scaled", batch_key="batch_str")
model = cytovi.CYTOVI(
    adata=adata_to_correct, n_hidden=par["n_hidden"], n_layers=par["n_layers"]
)
model.train()

# get batch corrected data
corrected_data = model.get_normalized_expression()

# have to add in the uncorrected markers as well
uncorrected_data = adata[:, markers_not_correct].layers["preprocessed"]

out_matrix = np.concatenate([corrected_data.to_numpy(), uncorrected_data], axis=1)
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

# run umap for quick check
# import scanpy as sc
# test_adata = out_adata.copy()
# test_adata = test_adata[:, markers_to_correct]
# test_adata.X = test_adata.layers["integrated"]
# test_adata.obs = adata.obs
# sc.pp.neighbors(test_adata, use_rep="X")
# sc.tl.umap(test_adata)
# sc.pl.umap(test_adata, color="batch")

print("Write output AnnData to file", flush=True)

out_adata.write_h5ad(par["output"], compression="gzip")
