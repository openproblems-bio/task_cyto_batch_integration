import sys

import anndata as ad

## VIASH START
par = {
    "input_unintegrated": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/unintegrated.h5ad",
    "output_integrated_split1": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/control_integrated_split1.h5ad",
    "output_integrated_split2": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/control_integrated_split2.h5ad",
}
meta = {"name": "shuffle_integration_by_cell_type", "resources_dir": "src/control_methods"}
## VIASH END

print("Importing helper functions", flush=True)
sys.path.append(meta["resources_dir"])
from utils import _randomize_features

print("Reading and preparing input files", flush=True)
adata = ad.read_h5ad(par["input_unintegrated"])
adata_split1 = adata[(adata.obs.is_control > 0) | (adata.obs.split == 1)].copy()
adata_split2 = adata[(adata.obs.is_control > 0) | (adata.obs.split == 2)].copy()

print("Randomise features - split 1", flush=True)
adata_split1.obs["batch_str"] = adata_split1.obs["batch"].astype(str)
integrated = _randomize_features(
    adata_split1.layers["preprocessed"]
)

# create new anndata
output_split1 = ad.AnnData(
    obs=adata_split1.obs[[]],
    var=adata_split1.var[[]],
    layers={"integrated": integrated},
    uns={
        "dataset_id": adata_split1.uns["dataset_id"],
        "method_id": meta["name"],
        "parameters": {},
    },
)

print("Randomise features - split 2", flush=True)
adata_split2.obs["batch_str"] = adata_split2.obs["batch"].astype(str)
integrated = _randomize_features(
    adata_split2.layers["preprocessed"]
)
# create new anndata
output_split2 = ad.AnnData(
    obs=adata_split2.obs[[]],
    var=adata_split2.var[[]],
    layers={"integrated": integrated},
    uns={
        "dataset_id": adata_split2.uns["dataset_id"],
        "method_id": meta["name"],
        "parameters": {},
    },
)

print("Write output AnnData to file", flush=True)
output_split1.write_h5ad(par["output_integrated_split1"], compression="gzip")
output_split2.write_h5ad(par["output_integrated_split2"], compression="gzip")