import sys

import anndata as ad

## VIASH START
par = {
    "input_unintegrated": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/unintegrated.h5ad",
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

print("Randomise features", flush=True)
integrated = _randomize_features(
    adata.layers["preprocessed"],
    partition=adata.obs["cell_type"],
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

all(x==y for x,y in zip(output.var_names, adata.var_names))

print("Write output AnnData to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
