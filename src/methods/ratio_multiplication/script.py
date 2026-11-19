import sys

import anndata as ad

## VIASH START
par = {
    "input": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/censored_split1.h5ad",
    "output": "output.h5ad",
    "percentile": 99.0,
}
meta = {"name": "ratio_multiplication", "resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from helper_functions import subset_nocontrols

print("Reading input files", flush=True)
adata = subset_nocontrols(ad.read_h5ad(par["input"]))

data = adata.to_df(layer="preprocessed").astype(float)
markers = adata.var_names[adata.var["to_correct"]]
in_batch1 = adata.obs["batch"].astype(str) == "1"
in_batch2 = adata.obs["batch"].astype(str) == "2"
q = par["percentile"] / 100

print("Correcting batch 2 to batch 1 via ratio multiplication", flush=True)
p1 = data.loc[in_batch1, markers].quantile(q)
p2 = data.loc[in_batch2, markers].quantile(q)
data.loc[in_batch2, markers] *= p1 / p2

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    layers={"integrated": data.to_numpy()},
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "method_id": meta["name"],
        "parameters": {"percentile": par["percentile"]},
    },
)
output.write_h5ad(par["output"], compression="gzip")
