import time

import anndata as ad
import numpy as np
import scvi
import torch
from scvi.external import cytovi
from sklearn.preprocessing import MinMaxScaler

## VIASH START
par = {
    "input": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_split2.h5ad",
    "output": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/output_cytovi_split2.h5ad",
    "n_hidden": 128,
    "n_layers": 1,
    "max_epochs": 1000,
    "train_size": 0.9,
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

print("Scaling data", flush=True)

# scale data. this will add a layer "scaled" to the anndata
cytovi.scale(
    adata=adata_to_correct,
    transformed_layer_key="preprocessed",
    batch_key="batch_str",
    scaled_layer_key="scaled",
    inplace=True,
    method="minmax",
)

# create minmax scaler object for one batch to use later for untransforming batch corrected data
print("Creating minmax scaler for untransforming data", flush=True)

# the following are taken from cytovi source code
feature_range = (0.0, 1.0)
feat_eps = 1e-6
feature_range = (feature_range[0] + feat_eps, feature_range[1] - feat_eps)

# get data from batch one
batch_one = (
    adata_to_correct[adata_to_correct.obs["batch_str"] == "1"]
    .layers["preprocessed"]
    .copy()
)
batch_one_scaler = MinMaxScaler(feature_range=feature_range)

print(
    "Fitting minmax scaler on batch one using feature range", feature_range, flush=True
)
batch_one_scaler.fit(batch_one)

print("Memory cleanup before training", flush=True)
del batch_one

print(
    f"Train CytoVI on {adata_to_correct.shape[0]} cells",
    flush=True,
)

cytovi.CYTOVI.setup_anndata(
    adata_to_correct,
    layer="scaled",
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
    max_epochs=par["max_epochs"],
    train_size=par["train_size"],
)
end = time.time()
print(f"Training took {end - start:.2f} seconds", flush=True)

# get batch corrected data
print("Calculating batch corrected data", flush=True)
corrected_data = model.get_normalized_expression(adata=adata_to_correct)

# have to save the columns to be able to reconstruct later
corrected_data_markers = corrected_data.columns.to_numpy()

print("Untransforming batch corrected data", flush=True)
# untransform data using batch one scaler
corrected_data = batch_one_scaler.inverse_transform(corrected_data.to_numpy())

# have to add in the uncorrected markers as well
uncorrected_data = adata[:, markers_not_correct].layers["preprocessed"].copy()

out_matrix = np.concatenate([corrected_data, uncorrected_data], axis=1)
out_var_idx = np.concatenate([corrected_data_markers, markers_not_correct])

# create new anndata
out_adata = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var.loc[out_var_idx][[]],
    layers={"integrated": out_matrix},
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "method_id": meta["name"],
        "parameters": par,
    },
)

# reorder var to match input
out_adata = out_adata[:, adata.var_names]

print("Write output AnnData to file", flush=True)

out_adata.write_h5ad(par["output"], compression="gzip")
