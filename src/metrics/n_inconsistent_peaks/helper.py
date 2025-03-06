import numpy as np
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

def get_kde_density(expression_array):
    ''' 
    Returns the density of the array using a gaussian kernel density estimation.
    
    Inputs:
    expression_array: array of values (1D) representing the marker expression

    Outputs:
    density: array of values representing the density of marker expression
    '''

    min_val = expression_array.min()
    max_val = expression_array.max()
    marker_values = np.reshape(expression_array, (1,-1)) # Reshape array for KDE
    kde = gaussian_kde(marker_values, bw_method=0.3)
    x_grid = np.linspace(min_val, max_val, 100)
    density = kde(x_grid)
    #Plot, for debugging
    # sns.scatterplot(x=x_grid, y=density)
    # plt.show()
    return density

def call_peaks(density):
    ''' 
    Returns the peaks of the density using scipy.signal.find_peaks.
    
    Inputs:
    density: array of values representing the density of marker expression

    Outputs:
    peaks: array of values representing the peaks of the density
    '''

    height_trsh = 0.05*density.max()
    prom_trsh = 0.1
                
    peaks, _ = find_peaks(density, 
                        prominence=prom_trsh,
                        height= height_trsh
                        )
    num_peaks = len(peaks)
        
    return num_peaks