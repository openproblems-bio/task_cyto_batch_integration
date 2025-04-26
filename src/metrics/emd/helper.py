import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance


def compute_emd(
    integrated_ct: ad.AnnData, 
    validation_ct: ad.AnnData, 
    markers_to_assess: list
) -> pd.DataFrame:
    """
    Calculate EMD metric

    Args:
        integrated_ct (ad.AnnData): batch integrated data
        validation_ct (ad.AnnData): validation data
        markers_to_assess (list): list of markers to compute EMD for

    Returns:
        pd.DataFrame: 1 row data frame where a column is a marker. Value is the EMD.
    """
    

    emd_vals = {}
        
    for marker in markers_to_assess:
        # marker = markers_to_assess[0]
        
        mexp_integrated = np.array(integrated_ct[:, marker].layers["integrated"]).flatten()
        mexp_validation = np.array(validation_ct[:, marker].layers["preprocessed"]).flatten()
        
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
    bin_edges = np.arange(-100, 100.1, 0.1)+0.0000001

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