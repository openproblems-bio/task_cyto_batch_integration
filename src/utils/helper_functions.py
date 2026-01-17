import anndata as ad
import numpy as np
import pandas as pd


def get_obs_var_for_integrated(
    s1_adata: ad.AnnData, s2_adata: ad.AnnData, u_adata: ad.AnnData
):
    """
    Fetch annotations (.var and .obs) from the unintegrated dataset to the integrated datasets (left and right).
    In the case of the control method 'perfect_integration', annotations are fetched only from batch 1.

    Inputs:
    s1_adata: AnnData object, split == 1 dataset
    s2_adata: AnnData object, split == 2 dataset
    u_adata: AnnData object, unintegrated dataset

    Outputs:
    s1_adata: AnnData object, split == 1 dataset with annotations
    s2_adata: AnnData object, split == 2 dataset with annotations
    """

    s1_adata.obs = u_adata.obs.loc[s1_adata.obs_names]
    s1_adata.var = u_adata.var.loc[s1_adata.var_names]
    s2_adata.obs = u_adata.obs.loc[s2_adata.obs_names]
    s2_adata.var = u_adata.var.loc[s2_adata.var_names]

    if s1_adata.uns["method_id"] == "perfect_integration":
        print(
            "Control method 'perfect_integration' detected. Changing batch labels for non-control cells in each split as they are all from batch 1."
        )

        # the idea is, we will change the batch labels in each split to reflect the correct donor/batch mapping
        # in unintegrated data.
        # e.g., if donor A is meant to come from batch 2 in split 1, we change the batch label of donor A in split 1 to 2.

        # Apply mapping to all non-control cells of split 1 and split 2 (+ change the split label)

        print("Old split 1 donor/batch mapping:", flush=True)
        print(
            s1_adata.obs.loc[:, ["donor", "batch", "is_control"]].value_counts(),
            flush=True,
        )

        split_dict_s1 = get_donor_batch_map(u_adata, split_of_interest=1)
        mask_s1 = s1_adata.obs.is_control == 0
        s1_adata.obs.loc[mask_s1, "batch"] = s1_adata.obs.loc[mask_s1, "donor"].map(
            split_dict_s1
        )
        s1_adata.obs.loc[mask_s1, "split"] = 1

        print("New split 1 donor/batch mapping:", flush=True)
        print(
            s1_adata.obs.loc[:, ["donor", "batch", "is_control"]].value_counts(),
            flush=True,
        )

        print("Old split 2 donor/batch mapping:", flush=True)
        print(
            s2_adata.obs.loc[:, ["donor", "batch", "is_control"]].value_counts(),
            flush=True,
        )

        split_dict_s2 = get_donor_batch_map(u_adata, split_of_interest=2)
        mask_s2 = s2_adata.obs.is_control == 0
        s2_adata.obs.loc[mask_s2, "batch"] = s2_adata.obs.loc[mask_s2, "donor"].map(
            split_dict_s2
        )
        s2_adata.obs.loc[mask_s2, "split"] = 2

        print("New split 2 donor/batch mapping:", flush=True)
        print(
            s2_adata.obs.loc[:, ["donor", "batch", "is_control"]].value_counts(),
            flush=True,
        )

    return s1_adata, s2_adata


def get_donor_batch_map(u_adata: ad.AnnData, split_of_interest: int) -> dict:
    """
    Create a dictionary that represent the correct donor/batch mapping for a split of interest.
    Note: This helper function is only meant to be used for 'perfect_integration' control method.

    Inputs:
    u_adata: AnnData object, unintegrated dataset
    split_of_interest: int, the split number from which the mapping is created

    Outputs:
    split_dict: Dictionary mapping donor/batch mapping for a split of interest
    """

    split_dict = (
        u_adata[
            (u_adata.obs.is_control == 0) & (u_adata.obs.split == split_of_interest)
        ]
        .obs.groupby("donor", observed=True)["batch"]
        .apply(pd.Series.unique)
        .apply(list)
        .to_dict()
    )

    # Safeguard to ensure that each donor has a unique batch in the split
    assert False not in [
        True if np.unique(el).size == 1 else False for el in split_dict.values()
    ], (
        "Donor/Batch mapping is ambiguous. Some donors have multiple batches in the same split."
    )

    # if batch column is string, return the batch values as string
    # else return as integer or whatever data type that matches.
    # TODO this is for backward compatibility, can be removed in future once we have
    # made sure all datasets have stringas batch labels.

    # this check the batch value (that is currently a list, but we just get one)
    # and see if it is string or not
    first_batch_val = list(split_dict.values())[0][0]
    if isinstance(first_batch_val, str):
        # transform the elements of the dictionary in strings
        split_dict = {k: str(v[0]) for k, v in split_dict.items()}
    elif isinstance(first_batch_val, (int, np.integer)):
        # only do int as otherwise it will be painful..
        # transform the elements of the dictionary in integers
        split_dict = {k: int(v[0]) for k, v in split_dict.items()}
    else:
        # error out as we don't know how to handle other data types
        raise ValueError(
            f"Batch value data type {type(first_batch_val)} not recognized. Only str and int are supported."
        )

    return split_dict


def subset_nocontrols(adata) -> ad.AnnData:
    """
    Subsets the anndata object to remove the control cells.
    Control cells are all the entries in adata.obs where is_control != 0.

    Inputs:
    adata: AnnData object

    Outputs:
    adata: AnnData object with cells from control samples removed
    """

    assert "is_control" in adata.obs.columns, (
        "The column 'is_control' is not present in the adata object."
    )

    # subset the adata to remove cells which is_control != 0
    adata = adata[adata.obs["is_control"] == 0].copy()

    return adata


def subset_markers_tocorrect(adata) -> ad.AnnData:
    """
    Subsets the anndata object to only include markers that need to be (or have been) corrected.
    This markers are all the entries in adata.var where to_correct == True.

    Inputs:
    adata: AnnData object

    Outputs:
    adata: AnnData object with only the markers to correct
    """

    adata = adata[:, adata.var["to_correct"]].copy()

    return adata


def remove_unlabelled(adata) -> ad.AnnData:
    """
    Subsets the anndata object to remove all cells where the marker is not labelled.
    This is determined by the column 'cell_type' in adata.obs.
    Particularly usefull when dealing with cell type specific metrics

    Inputs:
    adata: AnnData object

    Outputs:
    adata: AnnData object with only the labeled cells
    """

    is_unlabelled = adata.obs["cell_type"].str.lower().isin(["unlabelled", "unlabeled"])

    return adata[~is_unlabelled].copy()
