import itertools

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance

KEY_MEAN_EMD_GLOBAL = "mean_emd_global"
KEY_MAX_EMD_GLOBAL = "max_emd_global"
KEY_MEAN_EMD_CT = "mean_emd_ct"
KEY_MAX_EMD_CT = "max_emd_ct"
KEY_EMD_VERT_MAT_split1 = "emd_vert_mat_split1"
KEY_EMD_VERT_MAT_split2 = "emd_vert_mat_split2"
KEY_EMD_HORZ_PER_DONOR = "emd_horz_per_donor"


def calculate_vertical_emd(
    i_split1_adata: ad.AnnData, i_split2_adata: ad.AnnData, markers_to_assess: list
):
    """
    Compute vertical emd across every possible sample combination from the same biological group.

    Args:
        i_split1_adata (ad.AnnData): Left integrated adata.
        i_split2_adata (ad.AnnData): Right integrated adata.
        markers_to_assess (list): list of markers to compute EMD for.

    Returns:
        dict: a dictionary containing the following elements.
            "mean_emd_global": np.float32: mean emd value computed from a flattened data frame containing
                mean emd computed for every marker across all pairing two samples from the same group.
            "max_emd_global": np.float32: max emd value computed from a flattened data frame containing
                max emd computed for every marker across all pairing two samples from the same group.
            "mean_emd_ct": np.float32: mean emd value computed from a flattened data frame containing
                mean emd computed for every marker and cell type across all pairing two samples from the same group.
            "max_emd_ct": np.float32: max emd value computed from a flattened data frame containing
                max emd computed for every marker and cell type across all pairing two samples from the same group.
            "emd_wide_dfs": dict: each key is a marker. A value is yet another dictionary which value is
                a 2d matrix where each row/column is a pair of sample and a cell contains
                    emd value computed for the sample pair for either a given cell type
                    or global. The key for this dictionary then is either a given cell type
                    or global, depending on what the 2d matrix represents.
    """

    emd_split1_long, emd_split1_wide = get_vert_emd_for_integrated_adata(
        i_adata=i_split1_adata, markers_to_assess=markers_to_assess
    )

    emd_split2_long, emd_split2_wide = get_vert_emd_for_integrated_adata(
        i_adata=i_split2_adata, markers_to_assess=markers_to_assess
    )

    # safeguard
    mean_emd_global = np.nan
    max_emd_global = np.nan
    mean_emd_ct = np.nan
    max_emd_ct = np.nan

    # compute these only if we can.
    emd_long = []
    for df in [emd_split1_long, emd_split2_long]:
        if isinstance(df, pd.DataFrame):
            emd_long.append(df)

    if len(emd_long) > 0:
        emd_long = pd.concat(emd_long)

        # mean global emd across all sample combinations, markers, and splits
        mean_emd_global = np.nanmean(
            emd_long[emd_long["cell_type"] == "global"]
            .drop(columns=["cell_type", "first_sample", "second_sample"])
            .to_numpy()
            .flatten()
        )
        max_emd_global = np.nanmax(
            emd_long[emd_long["cell_type"] == "global"]
            .drop(columns=["cell_type", "first_sample", "second_sample"])
            .to_numpy()
            .flatten()
        )

        # mean cell type emd across all sample combinations, markers, and splits
        mean_emd_ct = np.nanmean(
            emd_long[emd_long["cell_type"] != "global"]
            .drop(columns=["cell_type", "first_sample", "second_sample"])
            .to_numpy()
            .flatten()
        )
        max_emd_ct = np.nanmax(
            emd_long[emd_long["cell_type"] != "global"]
            .drop(columns=["cell_type", "first_sample", "second_sample"])
            .to_numpy()
            .flatten()
        )

    return {
        KEY_MEAN_EMD_GLOBAL: mean_emd_global,
        KEY_MAX_EMD_GLOBAL: max_emd_global,
        KEY_MEAN_EMD_CT: mean_emd_ct,
        KEY_MAX_EMD_CT: max_emd_ct,
        KEY_EMD_VERT_MAT_split1: emd_split1_wide,
        KEY_EMD_VERT_MAT_split2: emd_split2_wide,
    }


def get_vert_emd_for_integrated_adata(i_adata: ad.AnnData, markers_to_assess: list):
    """
    Compute vertical emd across every possible sample combination from the same biological group.

    Args:
        i_split1_adata (ad.AnnData): An integrated adata.
        markers_to_assess (list): list of markers to compute EMD for.

    Returns:
        emd_dfs_ct: EMD per cell type,
        emd_dfs_global: EMD per donor,
        emd_wide_dfs: the break down of EMD per donor and cell type.
    """

    # calculate the global first, agnostic of cell type

    # comparing one sample against every other sample (one at a time).
    # get all samples for each group first
    sample_group_map = i_adata.obs.groupby("group", observed=True)["sample"].apply(
        lambda x: list(set(x))
    )
    # then get combinations of 2 for each group
    sample_combos = np.array(
        [list(itertools.combinations(x, 2)) for x in sample_group_map]
    ).reshape(-1, 2)

    if len(sample_combos) == 0:
        # this means the data processed do not have at least 2 samples per group.
        # thus it is impossible to compute this metric.

        print(
            f"{i_adata.uns['dataset_id']} from {i_adata.uns['method_id']} does not have"
            f" at least 2 samples per group. Skipping EMD vertical calculation."
        )

        return np.nan, np.nan

    cell_types = i_adata.obs["cell_type"].unique()

    emd_vals = []

    for sample_combo in sample_combos:
        # sample_combo = sample_combos[0]

        first_sample_adata = i_adata[i_adata.obs["sample"] == sample_combo[0]]
        second_sample_adata = i_adata[i_adata.obs["sample"] == sample_combo[1]]

        # global emd
        emd_df = compute_emd(
            left_sample=first_sample_adata,
            right_sample=second_sample_adata,
            markers_to_assess=markers_to_assess,
        )
        emd_df["cell_type"] = "global"
        emd_df["first_sample"] = sample_combo[0]
        emd_df["second_sample"] = sample_combo[1]
        emd_vals.append(emd_df)

        # emd per cell type
        for cell_type in cell_types:
            # cell_type = cell_types[0]
            first_sample_adata_ct = first_sample_adata[
                first_sample_adata.obs["cell_type"] == cell_type
            ]
            second_sample_adata_ct = second_sample_adata[
                second_sample_adata.obs["cell_type"] == cell_type
            ]

            # Do not calculate if we have less than 50 cells as it does not make sense.
            if first_sample_adata_ct.n_obs < 50 or second_sample_adata_ct.n_obs < 50:
                continue

            emd_df = compute_emd(
                left_sample=first_sample_adata_ct,
                right_sample=second_sample_adata_ct,
                markers_to_assess=markers_to_assess,
            )
            emd_df["cell_type"] = cell_type
            emd_df["first_sample"] = sample_combo[0]
            emd_df["second_sample"] = sample_combo[1]

            emd_vals.append(emd_df)

    # concatenate EMD values
    emd_vals = pd.concat(emd_vals)

    # prepare the data to draw the heatmap in cytonorm 2 supp paper.
    # 1 row/column = 1 sample, a cell is emd for a given marker
    # repeat for every marker assessed
    # note, only run this after calculating mean, otherwise you end up having to
    # remove the sample id columns.

    emd_wide = {}

    emd_types = emd_vals["cell_type"].unique()

    for marker in markers_to_assess:
        # marker = markers_to_assess[0]

        # remove unparsable characters like "/"
        marker_name = marker.replace("/", "_")
        # have to initialise the dictionary..
        emd_wide[marker_name] = {}

        for emd_type in emd_types:
            # ct = cell_types[0]
            emd_df = emd_vals[emd_vals["cell_type"] == emd_type]

            if emd_df.shape[0] > 0:
                # safeguard. Only pivot if we computed the emd.
                # This is a safeguard in case there is a rare cell type which we don't have
                # any samples with at least 50 cells for.
                emd_wide[marker_name][emd_type] = emd_df.pivot(
                    index="second_sample", columns="first_sample", values=marker
                )

    return emd_vals, emd_wide


def calculate_horizontal_emd(
    i_split1_adata: ad.AnnData,
    i_split2_adata: ad.AnnData,
    markers_to_assess: list,
    donor_list: list,
):
    """
    Compute horizontal emd across every pair samples from a given donor.

    Args:
        i_split1_adata (ad.AnnData): Left integrated adata.
        i_split2_adata (ad.AnnData): Right integrated adata.
        markers_to_assess (list): list of markers to compute EMD for.
        donor_list (list): list of donors to compute EMD for.

    Returns:
        dict: a dictionary containing the following elements.
            "mean_emd_global": np.float32: mean emd value computed from a flattened data frame containing
                mean emd computed for every marker across all pairs of samples from a given donor.
            "max_emd_global": np.float32: max emd value computed from a flattened data frame containing
                max emd computed for every marker across all pairs of samples from a given donor.
            "mean_emd_ct": np.float32: mean emd value computed from a flattened data frame containing
                mean emd computed for every marker and cell type across all pairs of samples from a given donor.
            "max_emd_ct": np.float32: max emd value computed from a flattened data frame containing
                max emd computed for every marker and cell type across all pairs of samples from a given donor.
            "emd_wide_dfs": pd.Dataframe showing emd value for each pair of sample
                and marker at either cell type level or global.
    """

    emd_per_donor_per_ct = []
    # global means agnostic of cell type labels
    emd_per_donor_global = []

    for donor in donor_list:
        # donor = donor_list[0]
        i_split1_donor = i_split1_adata[i_split1_adata.obs["donor"] == donor]
        i_split2_donor = i_split2_adata[i_split2_adata.obs["donor"] == donor]

        # safe check
        cell_type_not_in_both = np.setxor1d(
            i_split1_donor.obs["cell_type"].unique(),
            i_split2_donor.obs["cell_type"].unique(),
        )
        if len(cell_type_not_in_both) > 1:
            print(
                f"In donor {donor}: some cell types are in left integrated output"
                f" but not in right integrated output.\n"
                f"Cell types missing: {''.join(cell_type_not_in_both)}]n"
                f"Computing cell type EMD using just cell types common in both."
            )

        # assuming each cell type is present in both validation and integrated
        cell_types = np.intersect1d(
            i_split1_donor.obs["cell_type"].unique(),
            i_split2_donor.obs["cell_type"].unique(),
        )

        for cell_type in cell_types:
            # cell_type = cell_types[0]
            i_split1_ct = i_split1_donor[i_split1_donor.obs["cell_type"] == cell_type]
            i_split2_ct = i_split2_donor[i_split2_donor.obs["cell_type"] == cell_type]

            # Do not calculate if we have less than 50 cells as it does not make sense.
            if i_split1_ct.n_obs < 50 or i_split2_ct.n_obs < 50:
                print(
                    f"There are less than 50 cells for either left or right integrated "
                    f"data for donor {donor} and cell type {cell_type}.\n"
                    f"Skipping calculating EMD for this donor and cell type."
                )
                continue

            emd_df = compute_emd(
                left_sample=i_split1_ct,
                right_sample=i_split2_ct,
                markers_to_assess=markers_to_assess,
            )
            emd_df["cell_type"] = cell_type
            emd_df["donor"] = donor

            emd_per_donor_per_ct.append(emd_df)

        # calculate EMD when combining all cell types as well.
        emd_df = compute_emd(
            left_sample=i_split1_donor,
            right_sample=i_split2_donor,
            markers_to_assess=markers_to_assess,
        )
        emd_df["cell_type"] = "global"
        emd_df["donor"] = donor

        emd_per_donor_global.append(emd_df)

    emd_per_donor_per_ct = pd.concat(emd_per_donor_per_ct)
    emd_per_donor_global = pd.concat(emd_per_donor_global)

    # compute the mean and max per ct and for global.
    mean_emd_ct = np.nanmean(
        emd_per_donor_per_ct.drop(columns=["cell_type", "donor"]).values
    )
    max_emd_ct = np.nanmax(
        emd_per_donor_per_ct.drop(columns=["cell_type", "donor"]).values
    )

    mean_emd_global = np.nanmean(
        emd_per_donor_global.drop(columns=["cell_type", "donor"]).values
    )
    max_emd_global = np.nanmax(
        emd_per_donor_global.drop(columns=["cell_type", "donor"]).values
    )

    # concatenate the global and cell type emd
    emd_per_donor = pd.concat([emd_per_donor_per_ct, emd_per_donor_global])

    emd_per_donor.columns = emd_per_donor.columns.str.replace("/", "_", regex=False)

    return {
        KEY_MEAN_EMD_GLOBAL: mean_emd_global,
        KEY_MAX_EMD_GLOBAL: max_emd_global,
        KEY_MEAN_EMD_CT: mean_emd_ct,
        KEY_MAX_EMD_CT: max_emd_ct,
        KEY_EMD_HORZ_PER_DONOR: emd_per_donor,
    }


def compute_emd(
    left_sample: ad.AnnData, right_sample: ad.AnnData, markers_to_assess: list
) -> pd.DataFrame:
    """
    Calculate EMD metric

    Args:
        first_sample (ad.AnnData): data for the first sample to compute EMD for.
            For both horizontal and vertical EMD, this can be the sample from one donor that has been batch corrected.
        second_sample (ad.AnnData): data for the first sample to compute EMD for.
            For horizontal EMD, this will be the sample from the SAME donor (as in first_sample)
            that was not batch corrected, i.e., the one used in validation data.
            For vertical EMD, this will be just another sample that has been batch corrected.
        markers_to_assess (list): list of markers to compute EMD for.

    Returns:
        pd.DataFrame: 1 row data frame where a column is a marker. Value is the EMD.
    """

    emd_vals = {}

    for marker in markers_to_assess:
        # marker = markers_to_assess[0]

        mexp_split1 = np.array(left_sample[:, marker].layers["integrated"]).flatten()
        mexp_split2 = np.array(right_sample[:, marker].layers["integrated"]).flatten()

        left_values, left_weights = bin_array(mexp_split1)
        right_values, right_weights = bin_array(mexp_split2)

        # left_values (and right_values) are the explicit support (set of all possible bin values)
        # of the probability distribution left_weights (and right_weights).
        emd = wasserstein_distance(
            left_values, right_values, left_weights, right_weights
        )
        emd_vals[marker] = [emd]

    emd_df = pd.DataFrame.from_dict(emd_vals)

    return emd_df


def bin_array(values):
    """
    Bin values into probability distribution.

    Args:
        values (list): values to bin

    Returns:
        list: Bin indices - centre of the bin.
        list: Probability distribution of the input values.
    """

    # 2000 bins, the 0.0000001 is to avoid the left edge being included in the bin
    # (Mainly impacting 0 values)
    # range is set to -100 to 100 with the assumption that the range of values for each marker
    # will not exceed this
    bin_edges = np.arange(-100, 100.1, 0.1) + 0.0000001

    # using histogram retains the physical meaning of distances between bins,
    # such that moving mass from bin [-5.0, -4.9) to [-4.9, -4.8) has lower cost than
    # moving it to [99.9, 100.0)
    counts_per_bin, _ = np.histogram(values, bins=bin_edges)

    # this converts distribution of absolute marker values to probability distribution.
    # it allows subsequent EMD comparison between datasets of different sizes (number of cells).
    bin_probabilities = counts_per_bin / np.sum(counts_per_bin)

    # if bin_edges = [0,1,2,3,4,5]
    # bin_edges[:-1] will give you [0,1,2,3,4]
    # bin_edges[1:] will give you [1,2,3,4,5]
    # sum will sum each element up and divide by 2 will give you the centre of the bin
    # so for the 1st bin, (0+1)/2 = 0.5
    bin_indices = (bin_edges[:-1] + bin_edges[1:]) / 2

    # the 1st return value is the bin indices
    return bin_indices, bin_probabilities
