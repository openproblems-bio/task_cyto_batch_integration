import sys

import anndata as ad

## VIASH START
par = {
    "input_unintegrated": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/unintegrated.h5ad",
    "output": "resources_test/task_cyto_batch_integration/mouse_spleen_flow_cytometry_subset/precomputed_metric_info.h5ad",
}
meta = {
    "name": "precompute_metrics_info",
    "resources_dir": "src/utils",
}
sys.path.append("src/metrics/functional_marker_preservation")
## VIASH END

# Viash copies helper.py (from functional_marker_preservation) and
# helper_functions.py into resources_dir.
sys.path.append(meta["resources_dir"])

# this is the functional marker preservation helper.
from helper import (
    compute_cohens_d,
    compute_mean_expression_per_sample,
    filter_to_sig_pairs,
    get_significant_pairs,
    infer_groups,
    run_wilcoxon_group_tests,
)
from helper_functions import (
    preprocess_adata_for_evaluation,
)

print("Reading input files", flush=True)
unintegrated = ad.read_h5ad(par["input_unintegrated"])

unintegrated = preprocess_adata_for_evaluation(unintegrated)
unintegrated = unintegrated[:, unintegrated.var["marker_type"] == "functional"].copy()

# ── Functional marker preservation ──────────────────────────────────────────
# Inferred once here and reused by the metric for every method, so "group A"
# and "group B" (and so the sign of every Cohen's d) mean the same thing in
# every batch and split.
group_a, group_b = infer_groups(unintegrated)
print(f"Group A: {group_a}, Group B: {group_b}", flush=True)

print("Computing mean expression per sample", flush=True)

batches = sorted(unintegrated.obs["batch"].unique())

mean_expr_unintegrated = {}
for batch in batches:
    print(f"Processing unintegrated batch {batch}", flush=True)
    mean_expr_unintegrated[f"batch_{batch}"] = compute_mean_expression_per_sample(
        unintegrated[unintegrated.obs["batch"] == batch],
        layer="preprocessed",
    )

# Compare the mean expression of each functional marker in each cell type
# between the two groups using a Wilcoxon test, separately per batch. A pair
# enters the baseline only if it is significant in both batches.
print("Wilcoxon — computing unintegrated baseline", flush=True)

wilcoxon_unintegrated = {}
sig_per_batch = {}
for batch in batches:
    print(f"Processing unintegrated batch {batch}", flush=True)

    batch_key = f"batch_{batch}"

    results = run_wilcoxon_group_tests(
        mean_expr_unintegrated[batch_key], group_a, group_b
    )
    sig_per_batch[batch_key] = get_significant_pairs(results)
    wilcoxon_unintegrated[batch_key] = results
    print(
        f"{batch_key}: {len(sig_per_batch[batch_key])} significant (marker, cell_type) pairs",
        flush=True,
    )

sig_unintegrated = (
    sig_per_batch[f"batch_{batches[0]}"] & sig_per_batch[f"batch_{batches[1]}"]
)
print(
    f"Number of cell type/marker pairs significant in both batches: {len(sig_unintegrated)}",
    flush=True,
)

# precompute cohen D for sig pairs
print("Computing Cohen's d", flush=True)
cohens_d_unintegrated = {}

print("Filtering mean expression to significant pairs", flush=True)
mean_expr_unintegrated_sig = {
    f"batch_{batch}": filter_to_sig_pairs(
        mean_expr_unintegrated[f"batch_{batch}"], sig_unintegrated
    )
    for batch in batches
}

for batch in batches:
    batch_key = f"batch_{batch}"
    print(f"Computing cohen's d for {batch}", flush=True)
    cohens_d_unintegrated[batch_key] = compute_cohens_d(
        mean_expr_unintegrated_sig[batch_key], group_a, group_b
    )
    print(
        f"{batch_key}: {len(cohens_d_unintegrated[batch_key])} (marker, cell_type) pairs",
        flush=True,
    )

# Compute average Cohen's d across the two batches.
# Inner join keeps only pairs present in both batches.
print("Computing mean Cohen's d across batches", flush=True)
d_b1 = cohens_d_unintegrated[f"batch_{batches[0]}"][
    ["marker", "cell_type", "cohens_d"]
].rename(columns={"cohens_d": "cohens_d_b1"})
d_b2 = cohens_d_unintegrated[f"batch_{batches[1]}"][
    ["marker", "cell_type", "cohens_d"]
].rename(columns={"cohens_d": "cohens_d_b2"})
d_unintegrated_mean = d_b1.merge(d_b2, on=["marker", "cell_type"])
d_unintegrated_mean["cohens_d_mean_unintegrated"] = d_unintegrated_mean[
    ["cohens_d_b1", "cohens_d_b2"]
].mean(axis=1)

# Write everything the metric needs into uns.
print("Write output AnnData to file", flush=True)
ad.settings.allow_write_nullable_strings = True
output = ad.AnnData()
output.uns["functional_marker_preservation"] = {
    "mean_expr_unintegrated": mean_expr_unintegrated,
    "sig_pairs": list(sig_unintegrated),
    "wilcoxon_unintegrated": wilcoxon_unintegrated,
    "cohens_d_unintegrated": cohens_d_unintegrated,
    "cohens_d_mean_unintegrated": d_unintegrated_mean,
    "group_a": group_a,
    "group_b": group_b,
}

output.write_h5ad(par["output"], compression="gzip")
