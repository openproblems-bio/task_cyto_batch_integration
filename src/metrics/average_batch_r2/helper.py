import numpy as np
import pandas as pd
import anndata as ad


def concat_paired_samples(adata_s1: ad.AnnData,
                          adata_s2: ad.AnnData) -> pd.DataFrame:
    '''
    Concatenate split 1 and split 2 datasets to a single dataframe. Columns are markers, except last column ('batch')

    Inputs:
    adata_s1: AnnData object, split 1 integrated dataset
    adata_s2: AnnData object, split 2 integrated dataset

    Returns:
    df: pd.DataFrame with the split 1 and split 2 datasets concatenated
    '''
    assert adata_s1.shape[1] == adata_s2.shape[1], "The number of markers in the split 1 and split 2 datasets do not match"

    df_s1 = adata_s1.to_df(layer='integrated')
    df_s1['batch'] = adata_s1.obs['batch']

    df_s2 = adata_s2.to_df(layer='integrated')
    df_s2['batch'] = adata_s2.obs['batch']

    df = pd.concat([df_s1, df_s2], axis=0)

    return df

def fit_r2(df: pd.DataFrame, markername: str) -> float:
    ''''
    Fit a linear regression model to the marker expression data and calculate the R^2 value.
    
    Inputs:
    df: pd.DataFrame from `concat_paired_samples`
    markername: str, name of the marker to fit the R^2 value for
    
    Outputs:
    r2: float, The Rˆ2 represents the proportion of the variance explained by the batch variable in the marker expression.
    '''
    from sklearn.linear_model import LinearRegression

    data = df[[markername, 'batch']]
    X = pd.get_dummies(data['batch'], drop_first=True)
    Y = data[markername]
    model = LinearRegression().fit(X,Y)
    r2 = model.score(X,Y)

    # #Uncomment for debugging
    # import seaborn as sns
    # import matplotlib.pyplot as plt
    # sns.scatterplot(x='batch', y=markername, data=data)
    # plt.show()
    # print(markername,"r2 =",r2)

    return r2


def batch_r2(adata_s1: ad.AnnData, adata_s2: ad.AnnData) -> tuple[list, list]:
    '''
    Calculate the batch R^2 metric given 2 paired samples (split 1 and split 2).
    For each marker, the function calculates the R^2 value between the marker expression and batch covariate.
    Note: since adata_s1 and adata_s2 are paired samples, they have to come from the same donor.

    Inputs:
    adata_s1: AnnData object, split 1 integrated dataset
    adata_s2: AnnData object, split 2 integrated dataset

    Outputs:
    markers_r2: list of floats, R^2 values for each marker
    markerlist: list of str, marker
    '''

    assert np.unique(adata_s1.obs[ 'donor']) == np.unique(adata_s2.obs[ 'donor']), "The donors in the split 1 and split 2 datasets do not match"

    df = concat_paired_samples(adata_s1, adata_s2)

    markers_r2 = []
    markerlist = []
    for marker in df.columns:
        if marker != 'batch':
            r2 = fit_r2(df, marker)
            markers_r2.append(r2)
            markerlist.append(marker)
            

    return markers_r2,markerlist