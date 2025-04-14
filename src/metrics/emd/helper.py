import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance


def compute_emd(integrated_ct, validation_ct, markers_to_assess):
    emd_vals = {}
        
    for marker in markers_to_assess:
        # marker = markers_to_assess[0]
        
        mexp_integrated = np.array(integrated_ct[:, marker].layers["integrated"]).flatten()
        mexp_validation = np.array(validation_ct[:, marker].layers["preprocessed"]).flatten()
        
        i_values, i_weights = bin_array(mexp_integrated)
        v_values, v_weights = bin_array(mexp_validation)
        
        emd = wasserstein_distance(i_values, v_values, i_weights, v_weights)
        emd_vals[marker] = [emd]
        
    emd_df = pd.DataFrame.from_dict(emd_vals)
    
    return emd_df


def bin_array(values):
    ''''
    Input:
    - values (array) : array of values

    Returns:
    > a tuple with two arrays: the first array contains the binning, the second array contains the bin weights used to compute the EMD in the 'compute_emds_fromdict_ct' function
    '''
    
    # 2000 bins, the 0.0000001 is to avoid the left edge being included in the bin 
    # (Mainly impacting 0 values)
    bins = np.arange(-100, 100.1, 0.1)+0.0000001
    counts, _ = np.histogram(values, bins=bins)
    
    return range(0,2000), counts/sum(counts)