import anndata as ad


def get_obs_var_for_integrated(
    i_adata: ad.AnnData, v_adata: ad.AnnData, u_adata: ad.AnnData
) -> ad.AnnData:
    """
    Adds annotations (.var and .obs) from the unintegrated dataset to the integrated dataset.
    In the case of the control method 'perfect_integration', the function will fetch annotations from the validation dataset instead.

    Inputs:
    i_adata: AnnData object, batch-integrated dataset
    v_adata: AnnData object, validation dataset
    u_adata: AnnData object, unintegrated dataset

    Outputs:
    i_adata: AnnData object with .var and .obs added
    """
    import warnings

    if i_adata.uns["method_id"] == "perfect_integration_horizontal":
        assert (
            i_adata.shape[0] == v_adata.shape[0]
        ), "The number of cells in the integrated (perfect_integration_horizontal) and validation datasets do not match"
        i_adata.obs = v_adata.obs.loc[i_adata.obs_names]
        i_adata.var = v_adata.var.loc[i_adata.var_names]

    elif i_adata.uns["method_id"] == "perfect_integration_vertical":

        obs_adata = ad.concat([v_adata, u_adata])
        # subset to just batch 1
        # purposely hard code this so if the batch used for perfect integration vertical
        # changes, this will have to be purposely changed.
        obs_adata = obs_adata[obs_adata.obs["batch"] == 1]

        assert (
            i_adata.shape[0] == obs_adata.shape[0]
        ), "The number of cells in the integrated (perfect_integration_vertical) and validation + unintegrated datasets do not match"
        i_adata.obs = obs_adata.obs.loc[i_adata.obs_names]
        i_adata.var = v_adata.var.loc[i_adata.var_names]

    else:
        assert (
            i_adata.shape[0] == u_adata.shape[0]
        ), "The number of cells in the integrated and unintegrated datasets do not match"
        if False in list(i_adata.obs.index == u_adata.obs.index):
            warnings.warn(
                "The cell ordering in the integrated and unintegrated datasets do not match"
            )

        i_adata.obs = u_adata.obs.loc[i_adata.obs_names]
        i_adata.var = u_adata.var.loc[i_adata.var_names]

    return i_adata


def subset_nocontrols(adata) -> ad.AnnData:
    """
    Subsets the anndata object to remove the control cells.
    Control cells are all the entries in adata.obs where is_control != 0.

    Inputs:
    adata: AnnData object

    Outputs:
    adata: AnnData object with cells from control samples removed
    """

    assert (
        "is_control" in adata.obs.columns
    ), "The column 'is_control' is not present in the adata object."

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

    adata = adata[adata.obs["cell_type"].str.lower() != "unlabelled"].copy()
    adata = adata[adata.obs["cell_type"].str.lower() != "unlabeled"].copy()

    return adata
