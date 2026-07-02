import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def infer_groups(adata):
    """
    Determine the two biological group labels (e.g. WT/KO) from adata.obs["group"].

    Sorting makes group_a/group_b deterministic across runs and datasets — calling
    this once on the unintegrated data and threading the result through every
    run_wilcoxon_group_tests()/compute_cohens_d() call (rather than letting each
    call infer groups from its own subset via .unique(), which reflects row
    order and isn't guaranteed to agree across batches/splits) keeps "group A"
    and "group B" meaning the same thing everywhere, so a positive Cohen's d
    always means the same direction no matter which batch, split, or method run
    it came from.

    Args:
        adata: AnnData with an obs column 'group' containing exactly 2 unique values.

    Returns:
        tuple: (group_a, group_b), sorted alphabetically.
    """
    groups = sorted(adata.obs["group"].unique())
    if len(groups) != 2:
        raise ValueError(f"Expected exactly 2 biological groups, found: {groups}")
    return groups[0], groups[1]


def compute_mean_expression_per_sample(adata, layer, min_cells=5):
    """
    Compute mean expression of each functional marker per (cell_type, sample) pair.

    Args:
        adata: AnnData with obs columns ['cell_type', 'sample', 'group', 'batch']
            and the specified layer.
        layer (str): Layer key to use for expression values.
        min_cells (int): Minimum number of cells required for a (cell_type, sample)
            pair to be included. Pairs below this threshold are dropped.

    Returns:
        pd.DataFrame: Long-format DataFrame with columns
            ['cell_type', 'sample', 'group', 'batch', 'marker', 'mean_expr'].
    """
    X = np.array(adata.layers[layer])
    expr_df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    expr_df = expr_df.join(adata.obs[["cell_type", "sample", "group", "batch"]])

    # Drop (cell_type, sample) pairs that have too few cells for a reliable mean.
    # observed=True ensures only combinations that actually exist in the data are
    # returned — without it, pandas expands all combinations of categorical values,
    # producing NaN rows for (cell_type, sample) pairs that don't co-occur.
    counts = (
        expr_df.groupby(["cell_type", "sample"], observed=True).size().rename("n_cells")
    )
    expr_df = expr_df.join(counts, on=["cell_type", "sample"])
    expr_df = expr_df[expr_df["n_cells"] >= min_cells].drop(columns="n_cells")

    mean_df = (
        expr_df.groupby(["cell_type", "sample", "group", "batch"], observed=True)[
            adata.var_names.tolist()
        ]
        .mean()
        .reset_index()
    )

    return mean_df.melt(
        id_vars=["cell_type", "sample", "group", "batch"],
        var_name="marker",
        value_name="mean_expr",
    )


def run_wilcoxon_group_tests(mean_expr_df, group_a, group_b):
    """
    For each (marker, cell_type) pair, run a two-sided Wilcoxon rank-sum test
    comparing the two biological groups (e.g. WT vs KO) on per-sample mean
    expressions. Raw p-values are returned without FDR correction, as the small
    number of samples per group (typically 2-3) makes the minimum achievable
    p-value ~0.1, rendering FDR adjustment counterproductive.

    group_a and group_b must be determined once by the caller (see
    infer_groups()) and passed to every call, rather than re-inferred per
    batch/split — inferring locally from each subset's own unique() order can
    silently assign "group A" and "group B" to opposite biological groups
    across different calls. Use get_significant_pairs() to threshold results.

    Args:
        mean_expr_df (pd.DataFrame): Output of compute_mean_expression_per_sample.
            Must contain columns ['marker', 'cell_type', 'group', 'mean_expr'].
        group_a (str): Label of the first biological group.
        group_b (str): Label of the second biological group.

    Returns:
        pd.DataFrame: Columns ['marker', 'cell_type', 'group_a', 'group_b',
            'u_statistic', 'pval']. Pairs with fewer than 2 samples in either
            group are excluded.
    """
    # group_a/group_b are trusted as given (see infer_groups()), so this isn't
    # re-deriving them — it's a guard against mean_expr_df itself being wrong,
    # e.g. a third group label from a data bug, or group_a/group_b having been
    # computed from the wrong column upstream. Fails loudly here rather than
    # silently dropping the offending rows further down in the groupby loop.
    unexpected_groups = set(mean_expr_df["group"].unique()) - {group_a, group_b}
    if unexpected_groups:
        raise ValueError(
            f"mean_expr_df contains group(s) other than {group_a!r}/{group_b!r}: "
            f"{unexpected_groups}"
        )

    results = []
    # grp is the subset of rows for one (marker, cell_type) combination.
    # vals_a / vals_b are the per-sample mean expressions for each biological group.
    for (marker, cell_type), grp in mean_expr_df.groupby(
        ["marker", "cell_type"], observed=True
    ):
        vals_a = grp[grp["group"] == group_a]["mean_expr"].values
        vals_b = grp[grp["group"] == group_b]["mean_expr"].values

        # Need at least 2 samples per group for a meaningful rank-sum test.
        if len(vals_a) < 2 or len(vals_b) < 2:
            continue

        result = mannwhitneyu(vals_a, vals_b, alternative="two-sided")
        results.append(
            {
                "marker": marker,
                "cell_type": cell_type,
                "u_statistic": result.statistic,
                "pval": result.pvalue,
            }
        )

    if not results:
        return pd.DataFrame(
            columns=["marker", "cell_type", "group_a", "group_b", "u_statistic", "pval"]
        )

    results_df = pd.DataFrame(results)
    results_df["group_a"] = group_a
    results_df["group_b"] = group_b

    return results_df


def _cohens_d(vals_a, vals_b):
    """
    Compute Cohen's d between two arrays of values.

    d = (mean_A - mean_B) / pooled_SD

    Pooled SD weights each group's variance by its degrees of freedom:
        pooled_SD = sqrt(((n_A-1)*var_A + (n_B-1)*var_B) / (n_A + n_B - 2))

    Returns NaN when pooled_SD is zero (all values identical across both groups).

    Args:
        vals_a (np.ndarray): Values for group A.
        vals_b (np.ndarray): Values for group B.

    Returns:
        float: Cohen's d, or NaN if undefined.
    """
    n_a, n_b = len(vals_a), len(vals_b)
    denom = n_a + n_b - 2
    # setting ddof to 1 to get unbiased sample variance using Bessel’s correction
    # (dividing by n-1)
    s_a, s_b = vals_a.var(ddof=1), vals_b.var(ddof=1)
    pooled_sd = np.sqrt(((n_a - 1) * s_a + (n_b - 1) * s_b) / denom)

    if pooled_sd == 0:
        return np.nan

    mean_a, mean_b = vals_a.mean(), vals_b.mean()
    return (mean_a - mean_b) / pooled_sd


def compute_cohens_d(mean_expr_df, group_a, group_b):
    """
    For each (marker, cell_type) pair, compute Cohen's d between the two
    biological groups on per-sample mean expressions.

    Cohen's d = (mean_A - mean_B) / pooled_SD

    where pooled_SD = sqrt(((n_A-1)*var_A + (n_B-1)*var_B) / (n_A + n_B - 2))

    A positive d means group A has higher mean expression than group B.
    Returns NaN for a pair when pooled_SD is zero (all samples have identical
    expression) or when either group has fewer than 2 samples.

    group_a and group_b must be determined once by the caller (see
    infer_groups()) and passed to every call, rather than re-inferred per
    batch/split. Cohen's d is signed, so if "group A" silently swapped
    between calls, the sign of d would flip along with it — corrupting any
    average or comparison taken across batches/splits.

    Args:
        mean_expr_df (pd.DataFrame): Output of compute_mean_expression_per_sample.
            Must contain columns ['marker', 'cell_type', 'group', 'mean_expr'].
        group_a (str): Label of the first biological group.
        group_b (str): Label of the second biological group.

    Returns:
        pd.DataFrame: Columns ['marker', 'cell_type', 'group_a', 'group_b',
            'cohens_d']. Pairs with fewer than 2 samples in either group are
            excluded.
    """
    # group_a/group_b are trusted as given (see infer_groups()), so this isn't
    # re-deriving them — it's a guard against mean_expr_df itself being wrong,
    # e.g. a third group label from a data bug, or group_a/group_b having been
    # computed from the wrong column upstream. Fails loudly here rather than
    # silently dropping the offending rows further down in the groupby loop.
    unexpected_groups = set(mean_expr_df["group"].unique()) - {group_a, group_b}
    if unexpected_groups:
        raise ValueError(
            f"mean_expr_df contains group(s) other than {group_a!r}/{group_b!r}: "
            f"{unexpected_groups}"
        )

    results = []
    # grp is the subset of rows for one (marker, cell_type) combination.
    # vals_a / vals_b are the per-sample mean expressions for each biological group.
    for (marker, cell_type), grp in mean_expr_df.groupby(
        ["marker", "cell_type"], observed=True
    ):
        vals_a = grp[grp["group"] == group_a]["mean_expr"].values
        vals_b = grp[grp["group"] == group_b]["mean_expr"].values

        if len(vals_a) < 2 or len(vals_b) < 2:
            continue

        results.append(
            {
                "marker": marker,
                "cell_type": cell_type,
                "cohens_d": _cohens_d(vals_a, vals_b),
            }
        )

    if not results:
        return pd.DataFrame(
            columns=["marker", "cell_type", "group_a", "group_b", "cohens_d"]
        )

    df = pd.DataFrame(results)
    df["group_a"] = group_a
    df["group_b"] = group_b

    return df


def filter_to_sig_pairs(mean_expr_df, sig_pairs):
    """
    Keep only rows whose (marker, cell_type) combination is in sig_pairs.

    Args:
        mean_expr_df (pd.DataFrame): Must contain columns ['marker', 'cell_type'].
        sig_pairs (set): Set of significant (marker, cell_type) tuples to keep.

    Returns:
        pd.DataFrame: Filtered mean_expr_df.
    """
    is_kept = [
        (marker, cell_type) in sig_pairs
        for marker, cell_type in zip(mean_expr_df["marker"], mean_expr_df["cell_type"])
    ]
    return mean_expr_df[is_kept]


def get_significant_pairs(test_results_df, pval_threshold=0.1):
    """
    Extract the set of (marker, cell_type) tuples whose raw p-value is at or
    below the threshold.

    The default threshold of 0.1 reflects the minimum achievable two-sided
    p-value with n=3 samples per group in a Wilcoxon rank-sum test.

    Args:
        test_results_df (pd.DataFrame): Output of run_wilcoxon_group_tests.
        pval_threshold (float): Raw p-value threshold for significance.

    Returns:
        set of (marker, cell_type) tuples.
    """
    sig = test_results_df[test_results_df["pval"] <= pval_threshold]
    return set(zip(sig["marker"], sig["cell_type"]))
