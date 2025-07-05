import itertools

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance

KEY_MEAN_EMD_GLOBAL = "mean_emd_global"
KEY_MAX_EMD_GLOBAL = "max_emd_global"
KEY_MEAN_EMD_CT = "mean_emd_ct"
KEY_MAX_EMD_CT = "max_emd_ct"
KEY_EMD_VERT_MAT = "emd_vert_mat"
KEY_EMD_HORZ_PER_DONOR = "emd_horz_per_donor"


def calculate_vertical_emd(input_integrated: ad.AnnData, markers_to_assess: list):
    """
    Compute vertical emd across every possible sample combination from the same biological group.

    Args:
        input_integrated (ad.AnnData): Integrated adata.
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

    # calculate the global first, agnostic of cell type

    # comparing one sample against every other sample (one at a time).
    # get all samples for each group first
    sample_group_map = input_integrated.obs.groupby("group", observed=True)[
        "sample"
    ].apply(lambda x: list(set(x)))
    # then get combinations of 2 for each group
    sample_combos = np.array(
        [list(itertools.combinations(x, 2)) for x in sample_group_map]
    ).reshape(-1, 2)

    if len(sample_combos) == 0:
        # this means the data processed do not have at least 2 samples per group.
        # thus it is impossible to compute this metric.
        return {
            KEY_MEAN_EMD_GLOBAL: np.nan,
            KEY_MAX_EMD_GLOBAL: np.nan,
            KEY_MEAN_EMD_CT: np.nan,
            KEY_MAX_EMD_CT: np.nan,
            KEY_EMD_VERT_MAT: np.nan,
        }

    cell_types = input_integrated.obs["cell_type"].unique()

    emd_dfs_global = []
    emd_dfs_ct = []

    for sample_combo in sample_combos:
        # sample_combo = sample_combos[0]

        first_sample_adata = input_integrated[
            input_integrated.obs["sample"] == sample_combo[0]
        ]
        second_sample_adata = input_integrated[
            input_integrated.obs["sample"] == sample_combo[1]
        ]

        # global emd
        emd_df = compute_emd(
            first_sample=first_sample_adata,
            second_sample=second_sample_adata,
            markers_to_assess=markers_to_assess,
            first_sample_layer_name="integrated",
            second_sample_layer_name="integrated",
        )
        emd_df["first_sample"] = sample_combo[0]
        emd_df["second_sample"] = sample_combo[1]
        emd_dfs_global.append(emd_df)

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
                first_sample=first_sample_adata_ct,
                second_sample=second_sample_adata_ct,
                markers_to_assess=markers_to_assess,
                first_sample_layer_name="integrated",
                second_sample_layer_name="integrated",
            )
            emd_df["cell_type"] = cell_type
            emd_df["first_sample"] = sample_combo[0]
            emd_df["second_sample"] = sample_combo[1]

            emd_dfs_ct.append(emd_df)

    # process the global emd
    emd_dfs_global = pd.concat(emd_dfs_global)
    # aggregate into one metric by flattening all the values in the data frame
    # into one giant array and take a mean
    mean_emd_global = np.nanmean(
        emd_dfs_global.drop(columns=["first_sample", "second_sample"])
        .to_numpy()
        .flatten()
    )
    max_emd_global = np.nanmax(
        emd_dfs_global.drop(columns=["first_sample", "second_sample"])
        .to_numpy()
        .flatten()
    )

    # process the cell type specific emd
    emd_dfs_ct = pd.concat(emd_dfs_ct)
    # aggregate into one metric by flattening all the values in the data frame
    # into one giant array and take a mean

    mean_emd_ct = np.nanmean(
        emd_dfs_ct.drop(columns=["cell_type", "first_sample", "second_sample"])
        .to_numpy()
        .flatten()
    )
    max_emd_ct = np.nanmax(
        emd_dfs_ct.drop(columns=["cell_type", "first_sample", "second_sample"])
        .to_numpy()
        .flatten()
    )

    # prepare the data to draw the heatmap in cytonorm 2 supp paper.
    # 1 row/column = 1 sample, a cell is emd for a given marker
    # repeat for every marker assessed
    # note, only run this after calculating mean, otherwise you end up having to
    # remove the sample id columns.

    emd_wide_dfs = {}
    for marker in markers_to_assess:
        # marker = markers_to_assess[0]

        # start with global
        emd_wide = emd_dfs_global.pivot(
            index="second_sample", columns="first_sample", values=marker
        )
        emd_wide_dfs[marker] = {"global": emd_wide}

        # then per cell type
        for ct in cell_types:
            # ct = cell_types[0]
            emd_df = emd_dfs_ct[emd_dfs_ct["cell_type"] == ct]

            if emd_df.shape[0] > 0:
                # safeguard. Only pivot if we compute the emd for the cell type.
                # This is a safeguard in case there is a rare cell type which we don't have
                # any samples with at least 50 cells for.
                emd_wide = emd_df.pivot(
                    index="second_sample", columns="first_sample", values=marker
                )
                emd_wide_dfs[marker][ct] = emd_wide

    return {
        KEY_MEAN_EMD_GLOBAL: mean_emd_global,
        KEY_MAX_EMD_GLOBAL: max_emd_global,
        KEY_MEAN_EMD_CT: mean_emd_ct,
        KEY_MAX_EMD_CT: max_emd_ct,
        KEY_EMD_VERT_MAT: emd_wide_dfs,
    }


def calculate_horizontal_emd(
    input_integrated: ad.AnnData,
    input_validation: ad.AnnData,
    markers_to_assess: list,
    donor_list: list,
):
    """
    Compute horizontal emd across every pair samples from a given donor.

    Args:
        input_integrated (ad.AnnData): Integrated adata.
        input_validation (ad.AnnData): Validation adata.
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
        integrated_view = input_integrated[input_integrated.obs["donor"] == donor]
        validation_view = input_validation[input_validation.obs["donor"] == donor]

        # assuming each cell type is present in both validation and integrated
        cell_types = validation_view.obs["cell_type"].unique()

        for cell_type in cell_types:
            # cell_type = cell_types[0]
            integrated_ct = integrated_view[
                integrated_view.obs["cell_type"] == cell_type
            ]
            validation_ct = validation_view[
                validation_view.obs["cell_type"] == cell_type
            ]

            # Do not calculate if we have less than 50 cells as it does not make sense.
            if integrated_ct.n_obs < 50 or validation_ct.n_obs < 50:
                continue

            emd_df = compute_emd(
                first_sample=integrated_ct,
                second_sample=validation_ct,
                markers_to_assess=markers_to_assess,
                first_sample_layer_name="integrated",
                second_sample_layer_name="preprocessed",
            )
            emd_df["cell_type"] = cell_type
            emd_df["donor"] = donor

            emd_per_donor_per_ct.append(emd_df)

        # calculate EMD when combining all cell types as well.
        emd_df = compute_emd(
            first_sample=integrated_view,
            second_sample=validation_view,
            markers_to_assess=markers_to_assess,
            first_sample_layer_name="integrated",
            second_sample_layer_name="preprocessed",
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

    return {
        KEY_MEAN_EMD_GLOBAL: mean_emd_global,
        KEY_MAX_EMD_GLOBAL: max_emd_global,
        KEY_MEAN_EMD_CT: mean_emd_ct,
        KEY_MAX_EMD_CT: max_emd_ct,
        KEY_EMD_HORZ_PER_DONOR: emd_per_donor,
    }


def compute_emd(
    first_sample: ad.AnnData,
    second_sample: ad.AnnData,
    markers_to_assess: list,
    first_sample_layer_name: str,
    second_sample_layer_name: str,
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
        first_sample_layer_name (str): name of the layer to get the data out of.
            Should be integrated for horizontal EMD.
        second_sample_layer_name (str): name of the layer to get the data out of.
            Should be integrated or preprocessed for vertical or horizontal EMD respectively.

    Returns:
        pd.DataFrame: 1 row data frame where a column is a marker. Value is the EMD.
    """

    emd_vals = {}

    for marker in markers_to_assess:
        # marker = markers_to_assess[0]

        mexp_integrated = np.array(
            first_sample[:, marker].layers[first_sample_layer_name]
        ).flatten()
        mexp_validation = np.array(
            second_sample[:, marker].layers[second_sample_layer_name]
        ).flatten()

        i_values, i_weights = bin_array(mexp_integrated)
        v_values, v_weights = bin_array(mexp_validation)

        # i_values (and v_values) are the explicit support (set of all possible bin values)
        # of the probability distribution i_weights (and v_weights).
        emd = wasserstein_distance(i_values, v_values, i_weights, v_weights)
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
