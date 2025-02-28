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

    