import anndata as ad
import numpy as np
import scib_metrics as sm


def compute_ilisi_per_group(integrated: ad.AnnData, split_label: str):
    """
    Compute iLISI per biological group, skipping groups where batch is confounded.

    Args:
        integrated (ad.AnnData): Integrated AnnData with obs columns 'group' and 'batch'.
        split_label (str): Label for the split (used in print statements).

    Returns:
        ilisi_per_group (dict): Per-group normalised iLISI scalar values.
            Empty dict if all groups are skipped due to batch-group confounding.
        ilisi_per_cell_per_group (dict): Per-group arrays of raw per-cell iLISI values.
            Empty dict if all groups are skipped due to batch-group confounding.
    """
    groups = integrated.obs["group"].unique()
    ilisi_per_group = {}
    ilisi_per_cell_per_group = {}

    for group in groups:
        group_adata = integrated[integrated.obs["group"] == group]
        n_batches_in_group = len(group_adata.obs["batch"].unique())

        if n_batches_in_group < 2:
            print(
                f"{split_label} - Group {group}: batch is confounded by group "
                f"(only 1 batch present). Skipping.",
                flush=True,
            )
            continue

        knn = sm.nearest_neighbors.pynndescent(
            group_adata.layers["integrated"], n_neighbors=100, random_state=0
        )
        ilisi_per_cell = sm.lisi_knn(knn, group_adata.obs["batch"])
        ilisi = (np.nanmedian(ilisi_per_cell) - 1) / (n_batches_in_group - 1)

        ilisi_per_group[group] = ilisi
        ilisi_per_cell_per_group[group] = ilisi_per_cell

    if len(ilisi_per_group) == 0:
        print(
            f"{split_label}: batch is confounded by group in all groups. "
            f"Returning NaN for iLISI.",
            flush=True,
        )

    return ilisi_per_group, ilisi_per_cell_per_group


def compute_clisi(integrated: ad.AnnData):
    """
    Compute cLISI for an integrated AnnData object.

    Args:
        integrated (ad.AnnData): Integrated AnnData with obs column 'cell_type'.

    Returns:
        clisi (float): Normalised cLISI score.
        clisi_per_cell (np.ndarray): Raw per-cell cLISI values.
    """
    n_celltypes = len(integrated.obs["cell_type"].unique())
    knn = sm.nearest_neighbors.pynndescent(
        integrated.layers["integrated"], n_neighbors=100, random_state=0
    )
    clisi_per_cell = sm.lisi_knn(knn, integrated.obs["cell_type"])
    clisi = (n_celltypes - np.nanmedian(clisi_per_cell)) / (n_celltypes - 1)
    return clisi, clisi_per_cell
