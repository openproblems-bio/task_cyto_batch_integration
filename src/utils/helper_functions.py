import anndata as ad


def get_obs_for_integrated(
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
    
    adata = input_integrated.copy()
    
    # oh man i end up having to do this because otherwise some metrics which does 
    # pairwise comparison between samples from the same donor will give you error..
    if adata.uns['method_id'] == 'perfect_integration':
        
        adata.obs = input_validation.obs.loc[adata.obs_names]
        
        # i have to change the sample names as otherwise the metrics will get confused and
        # give NaN because it is expecting at least 2 samples per donor.
        sample_donor_map = input_unintegrated.obs.groupby(by='donor', observed=True)['sample'].unique().to_dict()
        
        adata.obs['sample'] = [sample_donor_map[x][0] for x in adata.obs['donor']]
    
    else:
        adata.obs = input_unintegrated.obs.loc[adata.obs_names]
        
    # re-arrange the var so validation and integrated have the same var order
    adata = adata[:, input_unintegrated.var_names]
    
    return(adata)

def get_obs_var_for_integrated(i_adata:ad.AnnData,
                               v_adata:ad.AnnData,
                               u_adata:ad.AnnData
                               ) -> ad.AnnData:
    '''
    Adds annotations (.var and .obs) from the unintegrated dataset to the integrated dataset.
    In the case of the control method 'perfect_integration', the function will fetch annotations from the validation dataset instead.
    
    Inputs:
    i_adata: AnnData object, batch-integrated dataset
    v_adata: AnnData object, validation dataset
    u_adata: AnnData object, unintegrated dataset
    
    Outputs:
    i_adata: AnnData object with .var and .obs added
    '''
    import warnings

    if i_adata.uns['method_id'] == 'perfect_integration':
        assert i_adata.shape[0] == v_adata.shape[0], "The number of cells in the integrated (perfect_integration) and validation datasets do not match"
        i_adata.obs = v_adata.obs.loc[i_adata.obs_names]
        i_adata.var = v_adata.var.loc[i_adata.var_names]

    else:
        assert i_adata.shape[0] == u_adata.shape[0], "The number of cells in the integrated and unintegrated datasets do not match"
        if False in list(i_adata.obs.index == u_adata.obs.index):
            warnings.warn("The cell ordering in the integrated and unintegrated datasets do not match")
        
        i_adata.obs = u_adata.obs.loc[i_adata.obs_names]
        i_adata.var = u_adata.var.loc[i_adata.var_names]

    return i_adata

def subset_nocontrols(adata)-> ad.AnnData:
    '''
    Subsets the anndata object to remove the control cells.
    Control cells are all the entries in adata.obs where is_control != 0.
    
    Inputs:
    adata: AnnData object

    Outputs:
    adata: AnnData object with cells from control samples removed 
    '''

    assert 'is_control' in adata.obs.columns, "The column 'is_control' is not present in the adata object."

    #subset the adata to remove cells which is_control != 0
    adata = adata[adata.obs['is_control'] == 0].copy()
    
    return adata

def subset_markers_tocorrect(adata)-> ad.AnnData:
    '''
    Subsets the anndata object to only include markers that need to be (or have been) corrected.
    This markers are all the entries in adata.var where to_correct == True.
    
    Inputs:
    adata: AnnData object

    Outputs:
    adata: AnnData object with only the markers to correct
    '''
    
    adata = adata[:,adata.var['to_correct']].copy()
    
    return adata
    