import anndata as ad
import numpy as np
from scvi.external import cytovi
from sklearn.cluster import KMeans

## VIASH START
par = {
    "input": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_split2.h5ad",
    "output": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/output_cytovi_split2.h5ad",
    "n_hidden": 128,
    "n_layers": 1,
    "n_clusters": 10,
    "subsample_fraction": 0.5,
}
meta = {"name": "cytovi"}
## VIASH END

print("Reading and preparing input files", flush=True)
adata = ad.read_h5ad(par["input"])

adata.obs["batch_str"] = adata.obs["batch"].astype(str)

markers_to_correct = adata.var[adata.var["to_correct"]].index.to_numpy()
markers_not_correct = adata.var[~adata.var["to_correct"]].index.to_numpy()

adata_to_correct = adata[:, markers_to_correct].copy()

print("Scaling data", flush=True)

# scale data. this will add a layer "scaled" to the anndata
cytovi.scale(
    adata=adata_to_correct,
    transformed_layer_key="preprocessed",
    batch_key="batch_str",
    inplace=True,
)

print("Clustering using k-means with k =", par["n_clusters"], flush=True)
# cluster data using Kmeans
adata_to_correct.obs["clusters"] = (
    KMeans(n_clusters=par["n_clusters"], random_state=0)
    .fit_predict(adata_to_correct.layers["scaled"])
    .astype(str)
)
# concatenate obs so we can use it for subsampling
adata_to_correct.obs["sample_cluster"] = (
    adata_to_correct.obs["sample"].astype(str) + "_" + adata_to_correct.obs["clusters"]
)
# subsample cells without replacement
print("Subsampling cells", flush=True)
subsampled_cells = adata_to_correct.obs.groupby("sample_cluster")[
    "sample_cluster"
].apply(lambda x: x.sample(n=round(len(x) * par["subsample_fraction"]), replace=False))
# need the cell id included in the subsample
subsampled_cells_idx = [x[1] for x in subsampled_cells.index.to_list()]

adata_subsampled = adata_to_correct[subsampled_cells_idx, :].copy()

print(
    f"Train CytoVI on subsampled data containing {adata_subsampled.shape[0]} cells",
    flush=True,
)

cytovi.CYTOVI.setup_anndata(adata_subsampled, layer="scaled", batch_key="batch_str")
model = cytovi.CYTOVI(
    adata=adata_subsampled, n_hidden=par["n_hidden"], n_layers=par["n_layers"]
)
model.train()

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
