import anndata as ad
import numpy as np
import scib_metrics as sm


def compute_ilisi(adata: ad.AnnData, n_batches: int):
    """
    Compute iLISI for an integrated AnnData object.

    Args:
        adata (ad.AnnData): Integrated AnnData with obs column 'batch' and
            layer 'integrated'.
        n_batches (int): Total number of batches across the dataset.

    Returns:
        ilisi (float): Normalised iLISI score in [0, 1].
            0 = cells are surrounded by only one batch (no mixing).
            1 = cells are perfectly mixed across all batches.
        ilisi_per_cell (np.ndarray): Raw per-cell LISI values before normalisation.
    """
    knn = sm.nearest_neighbors.pynndescent(
        adata.layers["integrated"], n_neighbors=100, random_state=42
    )
    ilisi_per_cell = sm.lisi_knn(knn, adata.obs["batch"])

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
