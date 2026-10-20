import numpy as np
import pandas as pd
from helper_functions import (
    remove_unlabelled,
    subset_markers_tocorrect,
    subset_nocontrols,
)
from ratio_consistent_peaks.helper import (
    call_peaks,
    get_kde_density,
    standardise_marker_expression,
)

# A (donor, cell type) is only scored if both splits have at least this many cells.
MIN_CELLS_PER_SPLIT = 100


def precompute(unintegrated):
    """
    Count the peaks of every marker in the unintegrated data, per donor and
    cell type, separately for split 1 and split 2.

    Uses the same filters as the ratio_consistent_peaks metric applies to the
    integrated splits, so the (donor, cell type, marker) combinations here are
    the ones the metric scores.

    Args:
        unintegrated: AnnData with the unintegrated data.

    Returns:
        dict with `unintegrated_peaks`: a DataFrame with one row per
        (donor, cell_type, marker) and columns n_peaks_s1, n_peaks_s2 (number
        of peaks in each split) and sd_s1, sd_s2 (standard deviation of the
        unscaled expression in each split, used by the metric to spot markers
        that lose all variance after integration).
    """
    unintegrated = subset_nocontrols(unintegrated)
    unintegrated = subset_markers_tocorrect(unintegrated)
    unintegrated = remove_unlabelled(unintegrated)

    unintegrated_s1 = unintegrated[unintegrated.obs["split"] == 1]
    unintegrated_s2 = unintegrated[unintegrated.obs["split"] == 2]

    peak_rows = []
    for donor in unintegrated_s1.obs["donor"].unique():
        print(f"Processing donor {donor}", flush=True)
        s1_donor = unintegrated_s1[unintegrated_s1.obs["donor"] == donor]
        s2_donor = unintegrated_s2[unintegrated_s2.obs["donor"] == donor]

        for celltype in s1_donor.obs["cell_type"].unique():
            s1_ct = s1_donor[s1_donor.obs["cell_type"] == celltype]
            s2_ct = s2_donor[s2_donor.obs["cell_type"] == celltype]

            if s1_ct.shape[0] < MIN_CELLS_PER_SPLIT or s2_ct.shape[0] < MIN_CELLS_PER_SPLIT:
                print(
                    f"Skipping celltype {celltype} and donor {donor}: "
                    f"{s1_ct.shape[0]} cells in split 1, {s2_ct.shape[0]} in split 2, "
                    f"need at least {MIN_CELLS_PER_SPLIT} in both.",
                    flush=True,
                )
                continue

            print(f"Processing celltype {celltype}", flush=True)
            for marker in s1_ct.var.index:
                # Standardise using the pooled mean and sd of both splits,
                # then count the peaks of each split's density.
                s1_unscaled = np.array(s1_ct[:, marker].layers["preprocessed"])
                s2_unscaled = np.array(s2_ct[:, marker].layers["preprocessed"])
                s1_scaled, s2_scaled = standardise_marker_expression(
                    s1_unscaled, s2_unscaled
                )

                peak_rows.append(
                    {
                        "donor": donor,
                        "cell_type": celltype,
                        "marker": marker,
                        "n_peaks_s1": call_peaks(get_kde_density(s1_scaled)),
                        "n_peaks_s2": call_peaks(get_kde_density(s2_scaled)),
                        "sd_s1": np.std(s1_unscaled),
                        "sd_s2": np.std(s2_unscaled),
                    }
                )

    unintegrated_peaks = pd.DataFrame(
        peak_rows,
        columns=[
            "donor",
            "cell_type",
            "marker",
            "n_peaks_s1",
            "n_peaks_s2",
            "sd_s1",
            "sd_s2",
        ],
    )
    return {"unintegrated_peaks": unintegrated_peaks}
