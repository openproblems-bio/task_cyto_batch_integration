import numpy as np
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

    # If the highest value is at the first bin, shift bins by one and adjust x_grid
    if np.argmax(density) == 0 and density.size > 1:
        print("Shifting KDE bins by one as the highest density is at the first bin.")
        # orig_x_grid = x_grid.copy()
        # recale the grid so we only have 99 bins and shift everything by one to the right..
        x_grid = np.linspace(min_val, max_val, 99)
        density = kde(x_grid)

        # Prepend a zero so the beginning, but remove the last value to keep size consistent
        # as otherwise we will end up with an extra bin...
        density = np.concatenate([[0.0], density])
        # Have to use actual grid spacing to keep uniform spacing in x_grid.
        # Can't just blindly add 1.
        step = (max_val - min_val) / (len(x_grid)) if len(x_grid) > 1 else 0.0
        x_grid = np.concatenate([[min_val - step], x_grid])

    if plot:
        # only needed when debugging locally, so don't import at module level
        import matplotlib.pyplot as plt
        import seaborn as sns

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
