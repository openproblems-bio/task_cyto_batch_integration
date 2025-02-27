import anndata as ad
import numpy as np


def concatenate_inputs(
    input_integrated: ad.AnnData,
    input_validation: ad.AnnData,
    input_unintegrated: ad.AnnData
) -> ad.AnnData:
    """
    Concatenate validation and integrated data.
    Unintegrated data is needed to embed the obs for batch integrated data.

    Args:
        input_integrated (ad.AnnData): batch integrated data
        input_validation (ad.AnnData): validation data
        input_unintegrated (ad.AnnData): unintegrated data

    Returns:
        ad.AnnData: concatenated data with actual marker expression stored in "data" layer
    """
    # add all obs columns in the unintegrated data to integrated data
    # loc should order the obs based on obs_names
    
    input_integrated.obs = input_unintegrated.obs.loc[input_integrated.obs_names]
    
    # re-arrange the var so validation and integrated have the same var order
    input_integrated = input_integrated[:, input_validation.var_names]

    input_concat = ad.concat([input_integrated, input_validation])
    input_concat.layers["data"] = np.concatenate(
    (
        input_integrated.layers["integrated"],
        input_validation.layers["preprocessed"],
    )
    )

    return input_concat