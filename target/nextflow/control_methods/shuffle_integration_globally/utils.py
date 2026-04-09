import numpy as np

def _randomize_features(X, partition=None):
    """
    Taken and adapted from opsca-v1:
    https://github.com/openproblems-bio/openproblems/blob/acf5c95a7306b819c4a13972783433d0a48f769b/openproblems/tasks/_batch_integration/_common/methods/baseline.py#L13
    """
    X_out = X.copy()
    if partition is None:
        partition = np.full(X.shape[0], 0)
    else:
        partition = np.asarray(partition)
    for partition_name in np.unique(partition):
        partition_idx = np.argwhere(partition == partition_name).flatten()
        X_out[partition_idx] = X[np.random.permutation(partition_idx)]
    return X_out
