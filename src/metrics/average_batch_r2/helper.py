import numpy as np
import pandas as pd
import anndata as ad


def concat_paired_samples(adata_i: ad.AnnData,
                          adata_v:ad.AnnData) -> pd.DataFrame:
    '''
    Concatenate the integrated and validation datasets to a single dataframe. Columns are markers, except last column ('batch')

    Inputs:
    adata_i: AnnData object, batch-integrated dataset
    adata_v: AnnData object, validation dataset

    Returns:
    df: pd.DataFrame with the integrated and validation datasets concatenated
    '''
    assert adata_i.shape[1] == adata_v.shape[1], "The number of markers in the integrated and validation datasets do not match"

    df_int = adata_i.to_df(layer='integrated')
    df_int['batch'] = adata_i.obs['batch']

    df_val = adata_v.to_df(layer='preprocessed')
    df_val['batch'] = adata_v.obs['batch']

    df = pd.concat([df_int, df_val], axis=0)

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


def batch_r2(adata_i: ad.AnnData,
                           adata_v:ad.AnnData) -> (list, list):
    '''
    Calculate the batch R^2 metric given 2 paired samples (integrated and validation).
    For each marker, the function calculates the R^2 value between the marker expression and batch covariate.
    Note: since adata_i and adata_v are paired samples, they have to come from the same donor.

    Inputs:
    adata_i: AnnData object, batch-integrated dataset
    adata_v: AnnData object, validation dataset

    Outputs:
    markers_r2: list of floats, R^2 values for each marker
    markerlist: list of str, marker
    '''
    
    assert np.unique(adata_i.obs[ 'donor']) == np.unique(adata_v.obs[ 'donor']), "The donors in the integrated and validation datasets do not match"

    df = concat_paired_samples(adata_i, adata_v)

    markers_r2 = []
    markerlist = []
    for marker in df.columns:
        if marker != 'batch':
            r2 = fit_r2(df, marker)
            markers_r2.append(r2)
            markerlist.append(marker)
            

    return markers_r2,markerlist