import sys

import anndata as ad
import numpy as np
import pandas as pd

## VIASH START
par = {'input_unintegrated': '/Users/putri.g/Documents/cytobenchmark/datasets/mouse_spleen_flow_cytometry/unintegrated.h5ad', 'input_integrated_split1': '/Users/putri.g/Documents/cytobenchmark/analysis/runs/run_2026-04-26_21-34-33/mouse_spleen_flow_cytometry/method_out/perfect_integration_split1.h5ad', 'input_integrated_split2': '/Users/putri.g/Documents/cytobenchmark/analysis/runs/run_2026-04-26_21-34-33/mouse_spleen_flow_cytometry/method_out/perfect_integration_split2.h5ad', 'output': '/Users/putri.g/Documents/cytobenchmark/analysis/data/fmp/2026-07-02/mouse_spleen_flow_cytometry/functional_marker_preservation_perfect_integration.h5ad'}
meta = {'name': 'functional_marker_preservation', 'resources_dir': '/Users/putri.g/Documents/GitHub/task_cyto_batch_integration/src/utils'}
sys.path.insert(0, '/Users/putri.g/Documents/GitHub/task_cyto_batch_integration/src/metrics/functional_marker_preservation')
## VIASH END

sys.path.append(meta["resources_dir"])

import helper as metric_helper
from helper_functions import (
    get_obs_var_for_integrated,
    preprocess_adata_for_evaluation,
)

ad.settings.allow_write_nullable_strings = True

print("Reading input files", flush=True)
unintegrated = ad.read_h5ad(par["input_unintegrated"])
integrated_s1 = ad.read_h5ad(par["input_integrated_split1"])
integrated_s2 = ad.read_h5ad(par["input_integrated_split2"])

print("Formatting input files", flush=True)
integrated_s1, integrated_s2 = get_obs_var_for_integrated(
    integrated_s1, integrated_s2, unintegrated
)

print("Filtering integrated and unintegrated data", flush=True)
# preprocess_adata_for_evaluation removes controls, unlabelled cells, and non-corrected
# markers. The additional functional marker filter is specific to this metric.
integrated_s1 = preprocess_adata_for_evaluation(integrated_s1)
integrated_s1 = integrated_s1[
    :, integrated_s1.var["marker_type"] == "functional"
].copy()

integrated_s2 = preprocess_adata_for_evaluation(integrated_s2)
integrated_s2 = integrated_s2[
    :, integrated_s2.var["marker_type"] == "functional"
].copy()

unintegrated = preprocess_adata_for_evaluation(unintegrated)
unintegrated = unintegrated[:, unintegrated.var["marker_type"] == "functional"].copy()

# Inferred once from the unintegrated data and reused everywhere below, so
# "group A" and "group B" — and therefore the sign of every Cohen's d value —
# mean the same thing regardless of which batch, split, or method run produced it.
group_a, group_b = metric_helper.infer_groups(unintegrated)
print(f"Group A: {group_a}, Group B: {group_b}", flush=True)

# ── Compute mean expression per sample ──────────────────────────────────────
# Done once per batch/split and reused for both Wilcoxon and Cohen's d.
print("Computing mean expression per sample", flush=True)

batches = unintegrated.obs["batch"].unique()

mean_expr_unintegrated = {}
for batch in batches:
    print(f"  Unintegrated batch {batch}", flush=True)
    mean_expr_unintegrated[batch] = metric_helper.compute_mean_expression_per_sample(
        unintegrated[unintegrated.obs["batch"] == batch], layer="preprocessed"
    )

print("  Integrated split 1", flush=True)
mean_expr_s1 = metric_helper.compute_mean_expression_per_sample(
    integrated_s1, layer="integrated"
)
print("  Integrated split 2", flush=True)
mean_expr_s2 = metric_helper.compute_mean_expression_per_sample(
    integrated_s2, layer="integrated"
)

# ── Wilcoxon metric ──────────────────────────────────────────────────────────
# Unintegrated baseline: group A vs B within each batch separately.
# A pair enters the baseline only if significant in BOTH batches (intersection),
# ensuring the group difference is not batch-specific.
print("Wilcoxon — computing unintegrated baseline", flush=True)

wilcoxon_unintegrated = {}
sig_per_batch = {}
for batch in batches:
    results = metric_helper.run_wilcoxon_group_tests(
        mean_expr_unintegrated[batch], group_a, group_b
    )
    sig_per_batch[batch] = metric_helper.get_significant_pairs(results)
    wilcoxon_unintegrated[str(batch)] = results.to_dict(orient="list")
    print(
        f"  Batch {batch}: {len(sig_per_batch[batch])} significant (marker, cell_type) pairs",
        flush=True,
    )

sig_unintegrated = sig_per_batch[batches[0]] & sig_per_batch[batches[1]]
print(f"  Significant in both batches: {len(sig_unintegrated)}", flush=True)

if len(sig_unintegrated) == 0:
    print(
        "WARNING: No (marker, cell_type) pairs are significant in both unintegrated batches. "
        "Wilcoxon metric will be NaN.",
        flush=True,
    )
    wilcoxon_metric = np.nan
    wilcoxon_s1 = {}
    wilcoxon_s2 = {}
    sig_s1 = set()
    sig_s2 = set()
    sig_integrated = set()
else:
    print("Wilcoxon — computing integrated splits", flush=True)
    results_s1 = metric_helper.run_wilcoxon_group_tests(mean_expr_s1, group_a, group_b)
    results_s2 = metric_helper.run_wilcoxon_group_tests(mean_expr_s2, group_a, group_b)

    # Only count pairs that were in the unintegrated baseline.
    sig_s1 = metric_helper.get_significant_pairs(results_s1) & sig_unintegrated
    sig_s2 = metric_helper.get_significant_pairs(results_s2) & sig_unintegrated
    wilcoxon_s1 = results_s1.to_dict(orient="list")
    wilcoxon_s2 = results_s2.to_dict(orient="list")

    # A pair is preserved only if significant in BOTH integrated splits.
    sig_integrated = sig_s1 & sig_s2
    wilcoxon_metric = len(sig_integrated) / len(sig_unintegrated)
    print(f"  Split 1 num of significant pairs: {len(sig_s1)}", flush=True)
    print(f"  Split 2 num of significant pairs: {len(sig_s2)}", flush=True)
    print(
        f"  Num of significant pairs (both splits): {len(sig_integrated)}", flush=True
    )

print(f"Wilcoxon metric: {wilcoxon_metric}", flush=True)

# ── Cohen's d metric ─────────────────────────────────────────────────────────
# Ground truth: Cohen's d computed per batch from unintegrated, then averaged across batches.
# Integrated: Cohen's d computed per split (batch-agnostic). Each split is compared
# individually to the ground truth — abs() is taken before averaging so that opposite-
# direction errors in split 1 and split 2 cannot cancel each other out.
# Metric: mean( mean(|d_s1 - d_unintegrated|, |d_s2 - d_unintegrated|) ) across all pairs.
#
# Only (marker, cell_type) pairs in sig_integrated — those the Wilcoxon test found
# significant in both integrated splits — are considered. Effect size drift is not
# meaningful for pairs where the group difference didn't survive integration.
print("Cohen's d — computing unintegrated ground truth", flush=True)

if len(sig_integrated) == 0:
    print(
        "WARNING: No (marker, cell_type) pairs are significant in both integrated "
        "splits. Cohen's d metric will be NaN.",
        flush=True,
    )
    cohens_d_metric = np.nan
    d_comparison = pd.DataFrame(
        columns=[
            "marker",
            "cell_type",
            "cohens_d_b1",
            "cohens_d_b2",
            "cohens_d_mean_unintegrated",
            "cohens_d_s1",
            "cohens_d_s2",
            "abs_diff_s1",
            "abs_diff_s2",
            "abs_mean",
        ]
    )
else:
    mean_expr_unintegrated_sig = {
        batch: metric_helper.filter_to_sig_pairs(mean_expr_unintegrated[batch], sig_integrated)
        for batch in batches
    }
    mean_expr_s1_sig = metric_helper.filter_to_sig_pairs(mean_expr_s1, sig_integrated)
    mean_expr_s2_sig = metric_helper.filter_to_sig_pairs(mean_expr_s2, sig_integrated)

    cohens_d_unintegrated = {}
    for batch in batches:
        cohens_d_unintegrated[batch] = metric_helper.compute_cohens_d(
            mean_expr_unintegrated_sig[batch], group_a, group_b
        )
        print(
            f"  Batch {batch}: {len(cohens_d_unintegrated[batch])} (marker, cell_type) pairs",
            flush=True,
        )

    # Average Cohen's d across the two batches.
    # Inner join keeps only pairs present in both batches.
    d_b1 = cohens_d_unintegrated[batches[0]][["marker", "cell_type", "cohens_d"]].rename(
        columns={"cohens_d": "cohens_d_b1"}
    )
    d_b2 = cohens_d_unintegrated[batches[1]][["marker", "cell_type", "cohens_d"]].rename(
        columns={"cohens_d": "cohens_d_b2"}
    )
    d_unintegrated_mean = d_b1.merge(d_b2, on=["marker", "cell_type"])
    d_unintegrated_mean["cohens_d_mean_unintegrated"] = d_unintegrated_mean[
        ["cohens_d_b1", "cohens_d_b2"]
    ].mean(axis=1)

    print("Cohen's d — computing integrated splits", flush=True)
    d_s1 = metric_helper.compute_cohens_d(mean_expr_s1_sig, group_a, group_b)
    d_s2 = metric_helper.compute_cohens_d(mean_expr_s2_sig, group_a, group_b)

    # Prepare per-split d tables for the comparison merge.
    d_split1 = d_s1[["marker", "cell_type", "cohens_d"]].rename(
        columns={"cohens_d": "cohens_d_s1"}
    )
    d_split2 = d_s2[["marker", "cell_type", "cohens_d"]].rename(
        columns={"cohens_d": "cohens_d_s2"}
    )
    d_s1_and_s2 = d_split1.merge(d_split2, on=["marker", "cell_type"])

    # Per pair: |d_s1 - d_unintegrated| and |d_s2 - d_unintegrated|, then average the two.
    d_comparison = d_unintegrated_mean.merge(
        d_s1_and_s2,
        on=["marker", "cell_type"],
    )
    d_comparison["abs_diff_s1"] = (
        d_comparison["cohens_d_s1"] - d_comparison["cohens_d_mean_unintegrated"]
    ).abs()
    d_comparison["abs_diff_s2"] = (
        d_comparison["cohens_d_s2"] - d_comparison["cohens_d_mean_unintegrated"]
    ).abs()
    d_comparison["abs_mean"] = d_comparison[["abs_diff_s1", "abs_diff_s2"]].mean(axis=1)

    cohens_d_metric = d_comparison["abs_mean"].mean()

print(f"Cohen's d metric (mean abs diff): {cohens_d_metric}", flush=True)

# ── Write output ─────────────────────────────────────────────────────────────
print("Write output AnnData to file", flush=True)
# uns layout:
#   metric_ids / metric_values: one entry per metric (wilcoxon, cohens_d).
#   group_a / group_b: the two biological groups compared everywhere in this
#     file (see infer_groups() in helper.py). A positive cohens_d means group_a
#     had higher mean expression than group_b for that (marker, cell_type) pair.
#
#   wilcoxon section:
#     sig_unintegrated / sig_split1 / sig_split2 / sig_integrated: sorted lists
#       of (marker, cell_type) tuples significant at each stage.
#     wilcoxon_results: raw test tables per comparison round with columns
#       [marker, cell_type, group_a, group_b, u_statistic, pval].
#
#   cohens_d section:
#     cohen_d_results: merged table with content of unintegrated and split1 and split2 d values,
#       restricted to (marker, cell_type) pairs in sig_integrated.
#     cohens_d_metric_all: mean of all abs_diff values across s1 and s2 for all pairs (not averaged per pair first).
#
uns = {
    "dataset_id": integrated_s1.uns["dataset_id"],
    "method_id": integrated_s1.uns["method_id"],
    "metric_ids": [
        "functional_marker_preservation_wilcoxon",
        "functional_marker_preservation_cohens_d",
    ],
    "metric_values": [wilcoxon_metric, cohens_d_metric],
    "group_a": group_a,
    "group_b": group_b,
    # Wilcoxon
    "sig_unintegrated": sorted(sig_unintegrated),
    "sig_split1": sorted(sig_s1),
    "sig_split2": sorted(sig_s2),
    "sig_integrated": sorted(sig_integrated),
    "wilcoxon_results": {
        "unintegrated": wilcoxon_unintegrated,
        "split_1": wilcoxon_s1,
        "split_2": wilcoxon_s2,
    },
    # Cohen's d
    "cohens_d_results": d_comparison.to_dict(orient="list"),
}

output = ad.AnnData(uns=uns)
output.write_h5ad(par["output"], compression="gzip")

print(f"Anndata written to {par['output']}", flush=True)
