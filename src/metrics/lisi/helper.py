import anndata as ad
import numpy as np
import scib_metrics as sm


def _check_batch_group_confounding(adata: ad.AnnData) -> bool:
    """
    Check whether batch is confounded by group across the dataset.

    A group is considered confounded when all of its cells come from a single
    batch, meaning LISI cannot detect cross-batch mixing within that group.
    Prints a warning listing all confounded groups if any are found.

    Args:
        adata (ad.AnnData): AnnData with obs columns 'group' and 'batch'.

    Returns:
        bool: True if any group is confounded (iLISI is undefined),
            False if all groups span multiple batches.
    """
    groups = adata.obs["group"].unique()
    confounded_groups = []

    for group in groups:
        n_batches_in_group = len(
            adata[adata.obs["group"] == group].obs["batch"].unique()
        )
        if n_batches_in_group < 2:
            confounded_groups.append(group)

    any_confounding_groups = len(confounded_groups) > 0

    if any_confounding_groups:
        print(
            f"Batch is confounded by group in the following groups: {confounded_groups}. "
            "Returning NaN for iLISI.",
            flush=True,
        )
    return any_confounding_groups


def compute_ilisi(adata: ad.AnnData, n_batches: int):
    """
    Compute iLISI for an integrated AnnData object.

    Skips computation and returns NaN if batch is confounded by group in any
    group (i.e., at least one group has cells from only one batch).

    Args:
        adata (ad.AnnData): Integrated AnnData with obs columns 'batch' and
            'group', and layer 'integrated'.
        n_batches (int): Total number of batches across the dataset.

    Returns:
        ilisi (float): Normalised iLISI score in [0, 1].
            0 = cells are surrounded by only one batch (no mixing).
            1 = cells are perfectly mixed across all batches.
            NaN if batch is confounded by group in any group.
        ilisi_per_cell (np.ndarray): Raw per-cell LISI values before normalisation.
            Empty array if computation was skipped.
    """
    if _check_batch_group_confounding(adata):
        return np.nan, np.array([])

    knn = sm.nearest_neighbors.pynndescent(
        adata.layers["integrated"], n_neighbors=100, random_state=42
    )
    ilisi_per_cell = sm.lisi_knn(knn, adata.obs["batch"])

    # Raw LISI ranges from 1 (all neighbours from one batch) to n_batches
    # (neighbours uniformly distributed across batches). Subtract 1 and divide
    # by (n_batches - 1) to normalise to [0, 1].
    ilisi = (np.nanmedian(ilisi_per_cell) - 1) / (n_batches - 1)

    return ilisi, ilisi_per_cell


def compute_clisi(integrated: ad.AnnData):
    """
    Compute cLISI for an integrated AnnData object.

    Args:
        integrated (ad.AnnData): Integrated AnnData with obs column 'cell_type'
            and layer 'integrated'.

    Returns:
        clisi (float): Normalised cLISI score in [0, 1].
            0 = cell types are fully mixed in the neighbourhood (bad).
            1 = each cell's neighbourhood contains only its own cell type (good).
        clisi_per_cell (np.ndarray): Raw per-cell LISI values before normalisation.
    """
    n_celltypes = len(integrated.obs["cell_type"].unique())
    knn = sm.nearest_neighbors.pynndescent(
        integrated.layers["integrated"], n_neighbors=100, random_state=0
    )
    clisi_per_cell = sm.lisi_knn(knn, integrated.obs["cell_type"])

    clisi = (n_celltypes - np.nanmedian(clisi_per_cell)) / (n_celltypes - 1)

    return clisi, clisi_per_cell
