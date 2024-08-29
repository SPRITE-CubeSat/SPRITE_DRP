"""
tools to reconstruct an SNR image from brightness profiles
extracted using SPRITE push-broom mapping
"""

#######################
# Third Party Imports #
#######################

# Numpy
import numpy as np

# Scipy
from scipy.interpolate import interp1d

#################
# Local Imports #
#################

from sprite_sim_config import *

global x_pix, y_pix, exposure_time, quantum_efficiency


#################################################
# Convert 1D brightness profiles into 2D arrays #
#     Photon values for each pixel of slit      #
#################################################

def convert_brightness_profile_to_2d_array(e_per_s):
    """
    Parameters:
    -----------
    e_per_s: numpy array
        electron per second count for a given coordinate
        
    Returns
    -------
    interpolated_counts_2d: numpy array
        SNR brightness values in shape of slit
        for a singular profile
    
    """
    # Convert from electrons/s to photons (assume quantum efficiency of 0.5)
    photons = (e_per_s * exposure_time) / quantum_efficiency

    # Remove background noise
    outliers = np.where(photons > (5 * np.median(photons)))[0]
    photons[outliers] = np.mean(photons)

    ###############
    # Interpolate #
    ###############

    # Reverse order
    photons = photons[::-1]

    # ORIGINAL LENGTH OF PROFILE COUNTS
    original_x = np.arange(len(photons))

    # Slit y - dimensions (profile extracted along y-axis)
    target_x = np.linspace(0, len(photons) - 1, y_pix)

    # Interpolate 1D data to be length of slit y - pixels
    interpolating_function = interp1d(original_x, photons, kind='linear')
    interpolated_counts_1d = interpolating_function(target_x)

    # Create a 2D array by repeating the 1D interpolated data 40 times (number of x-pixels per slit)
    interpolated_counts_2d = np.tile(interpolated_counts_1d, (x_pix, 1))

    return interpolated_counts_2d


##########################################
# Create stacked array of SNR image data #
##########################################

def stack_slit_data_to_reconstruct_supernova_image(brightness_arrays):
    """
    Parameters:
    -----------
    brightness_arrays:
        list of brightness data arrays for each slit
        in a push-broom mapping
        
    Returns
    -------
    flipped_array: numpy array
        SNR brightness values in shape of image

    """

    ##########################
    # combine profile arrays # 
    ##########################

    # If we called in brightness profiles
    # Vertically combine profiles to form one 2D array
    combined_array = np.vstack(brightness_arrays)

    # Flip array to match SNR image orientation
    flipped_array = combined_array.T

    return flipped_array
