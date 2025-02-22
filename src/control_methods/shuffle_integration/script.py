import anndata as ad
import sys

## VIASH START
par = {
    "input_unintegrated": "resources_test/task_cyto_batch_integration/starter_file/unintegrated_censored.h5ad",
    "output": "output.h5ad",
}
meta = {"name": "harmonypy"}
## VIASH END

print("Importing helper functions", flush=True)
sys.path.append(meta["resources_dir"])
from utils import _randomize_features

print("Reading and preparing input files", flush=True)
adata = ad.read_h5ad(par["input_unintegrated"])

adata.obs["batch_str"] = adata.obs["batch"].astype(str)

markers_to_correct = adata.var[adata.var["to_correct"]].index.to_numpy()

adata = adata[:, markers_to_correct]

print("Randomise features", flush=True)
integrated = _randomize_features(
    adata.layers["preprocessed"]
)

# create new anndata
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    layers={"integrated": integrated},
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "method_id": meta["name"],
        "parameters": {},
    },
)

print("Write output AnnData to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
