import sys

import anndata as ad

## VIASH START
par = {
    "input_unintegrated": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/unintegrated.h5ad",
    "output": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/precomputed_python_metric.h5ad",
}
meta = {
    "name": "precompute_python_metric",
    "resources_dir": "src/utils",
}
# the per-metric precompute modules import each metric's helper.py as
# <metric_name>.helper, so src/metrics needs to be importable.
sys.path.append("src/metrics")
sys.path.append("src/data_processors/precompute_python_metric")
## VIASH END

# Viash copies the per-metric precompute modules, helper_functions.py and each
# metric's helper.py (into a folder named after the metric) into resources_dir.
sys.path.append(meta["resources_dir"])

import precompute_functional_marker_preservation
import precompute_ratio_consistent_peaks

# One precompute module per metric. Each module filters the unintegrated data
# itself and returns a dict, stored in uns under the metric's name. To change
# what's precomputed for one metric, only its module needs editing.
precompute_per_metric = {
    "functional_marker_preservation": precompute_functional_marker_preservation.precompute,
    "ratio_consistent_peaks": precompute_ratio_consistent_peaks.precompute,
}

print("Reading input files", flush=True)
unintegrated = ad.read_h5ad(par["input_unintegrated"])

output = ad.AnnData()
for metric_name, precompute in precompute_per_metric.items():
    print(f"── Precomputing for {metric_name} ──", flush=True)
    output.uns[metric_name] = precompute(unintegrated)

print("Write output AnnData to file", flush=True)
ad.settings.allow_write_nullable_strings = True
output.write_h5ad(par["output"], compression="gzip")
