import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from ripser import ripser
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde


def standardise_marker_expression(dist_1, dist_2):
    """
    Standardises the marker expression values from two distributions.

    Inputs:
    dist_1: array of values (1D) representing the marker expression from distribution 1
    dist_2: array of values (1D) representing the marker expression from distribution 2

    Outputs:
    std_dist_1: array of standardised values for distribution 1
    std_dist_2: array of standardised values for distribution 2
    """

    pooled = np.concatenate([dist_1, dist_2])
    mu, sd = pooled.mean(), pooled.std()
    std_dist_1 = (dist_1 - mu) / (sd)
    std_dist_2 = (dist_2 - mu) / (sd)

    return std_dist_1, std_dist_2


def get_kde_density(expression_array, return_xgrid=False, plot=False):
    """
    Returns the density of the array using a gaussian kernel density estimation.

    Inputs:
    expression_array: array of values (1D) representing the marker expression
    return_xgrid: boolean, if True, also return the x_grid values used for density estimation
    plot: boolean, if True, plot the density estimation

    Outputs:
    density: array of values representing the density of marker expression
    x_grid (optional): array of x values where the density is evaluated
    """

    min_val = expression_array.min()
    max_val = expression_array.max()
    marker_values = np.reshape(expression_array, (1, -1))  # Reshape array for KDE
    kde = gaussian_kde(marker_values, bw_method="scott")
    x_grid = np.linspace(min_val, max_val, 100)
    density = kde(x_grid)

    if plot:
        fig, ax = plt.subplots()
        sns.scatterplot(x=x_grid, y=density, ax=ax)
        ax.set_title("KDE Density Estimation")
        ax.set_xlabel("Marker Expression")
        ax.set_ylabel("Density")
        fig.tight_layout()
        fig.show()

    if return_xgrid:
        # handy for plotting later on and maybe even save in the AnnData object
        return density, x_grid
    else:
        return density


def call_peaks(density):
    """
    Returns the peaks of the density using scipy.signal.find_peaks.

    Inputs:
    density: array of values representing the density of marker expression

    Outputs:
    peaks: array of values representing the peaks of the density
    """

    height_trsh = 0.1
    prom_trsh = 0.01

    peaks, _ = find_peaks(density, prominence=prom_trsh, height=height_trsh)
    num_peaks = len(peaks)

    return num_peaks


def persistent_peak_count(ys, persistence_cutoff=0.08):
    """
    Counts robust peaks in a 1D dataset using persistent homology.

    Args:
        ys (np.ndarray): KDE of a marker expression (1D array)
        persistence_cutoff (float): a threshold that decides which peaks are “significant enough” to count.
            A large persistence peak survives over many levels of smoothing (i.e. a strong, real peak).
            A small persistence peak quickly merges into a neighbor — likely noise.
            0.01: very low threshold counts even weak bumps as peaks
            0.05: moderate (default) counts clearly separated peaks
            0.1–0.2: high threshold counts only strong, dominant peaks
            Default to 0.08 to biased towards strong peaks but not overly.

    Returns:
        int: number of significant peaks
    """

    # Invert to turn peaks into "holes" for 0D persistence
    Y = -ys.reshape(-1, 1)
    diagram = ripser(Y, maxdim=0)["dgms"][0]
    persistence = diagram[:, 1] - diagram[:, 0]

    # Define significance threshold relative to data range
    threshold = persistence_cutoff * np.ptp(ys)
    n_peaks = np.sum(persistence > threshold)
    return n_peaks
