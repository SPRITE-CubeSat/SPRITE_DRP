"""
Tools to add shot noise and background to a 2D spectrum

** background sources include geo-coronal lyman alpha emission and dark noise **
** does not include geo-coronal noise for 1 SN target yet **
"""

###########
# imports #
###########

# Third-party Imports
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import scipy
from scipy.signal import gaussian

# Local Imports
from bin_and_convert_data import convert_flux_to_photons
from create_2d_supernova_spec import create_gaussian, convolve_data

############################
# Aperture/slit constants #
############################

from sprite_sim_config import *
global row_pixels, col_pixels, dark_rate, ang_range, exposure_time


##################
# Add shot noise #
##################
def generate_noisy_image(convolved_2d_spec, slit_data):
    
    """
    Takes the convolved 2D spectrum and adds random shot noise
    and background from dark rate and geo-coronal lyman alpha emission

    Parameters
    ----------
    convolved_2d_spec: 2D numpy array
        convolved and normalized
        SNR / star data

    slit_data: 2D numpy array
        slit data file with 760 rows and 1380 cols
        detector data

    Returns
    -------
    spec_2d_with_background: 2D numpy array
        contains SNR and star photon data with
        noise and background
    """
        
    # Go through each photon value in the 2D spec
    # Add noise randomly
    
    # 760 rows
    for row in range(row_pixels):
        # 1380 columns
        for col in range(col_pixels):
            # Make a noisy image
            photon_val = np.abs(convolved_2d_spec[row, col])
            
            convolved_2d_spec[row, col] = np.random.normal(photon_val, np.sqrt(photon_val))

    # New variable to differentiate between ideal 2D spectrum and noisy 2D spectrum
    noisy_2d_spectrum = convolved_2d_spec
     
    # Generate 2D spectrum with background noise
    spec_2d_with_background = add_background(noisy_2d_spectrum, slit_data)
    
    return spec_2d_with_background


############################
# generate dark background #
############################
def gen_dark_bg(noisy_2d_spectrum):
    """
    Adds dark counts to the 2D spectrum

    Parameters
    ----------
    noisy_2d_spectrum: 2D numpy array
        convolved and normalized SN data with
        added shot noise

    Returns
    -------
    dark_bg_spec: 2D array
        contains photon count data with background from
        dark rate added to N random pixels where
        count_num = darkRate * exposureTime
    """
    
    # Determine number of counts to add to image
    count_num = dark_rate * exposure_time
    
    # Initialize 2D array of photon counts to add dark counts to
    dark_bg = noisy_2d_spectrum

    # Randomly choose (count_num) x-positions to add 1 count to
    x_positions = np.random.randint(0, row_pixels-1, count_num)
    
    # Randomly choose (count_num) y-positions to add 1 count to
    y_positions = np.random.randint(0, col_pixels-1, count_num)
    
    # Iterate through x and y values
    for count in range(count_num):
        
        # Pixel positions
        x = x_positions[count]
        y = y_positions[count]
        
        # Add count to image
        dark_bg[x, y] += 1
        
    return dark_bg


###############################
# add geo-coronal lyA emission #
###############################
def gen_lyman_bg(slit_data):
    
    """
    Adds geo-coronal lyman alpha emission
    background to 2D slit

    Parameters
    ----------
    slit_data: 2D numpy array
        slit data file with 760 rows and 1380 cols
        of detector data

    Returns
    -------
    lyman_bg: 2D array
        contains photon count data with background from
        geo-coronal lyman alpha emission
    """

    # Wavelength to center background around
    lym_alpha_wave = 1217                     # Angstroms

    # Standard deviation for lyman cont. background
    lym_std_dev = 30                          # Pixels

    ###################################
    # determine background brightness #
    ###################################
    
    # Lyman background brightness
    lym_brightness = 0.0094                   # Counts/cm^2/arcsecond^2/s

    # Get SPRITE area at lyman wavelength      # Cm^2
    lym_area = 14.3303
    
    # Get brightness for entire slit in photon counts
    # Multiply by effective area, slit area of 1800 arcsecond^2, 100 arcsecond^2, exposure time
    lym_photons = lym_brightness * lym_area * 1800 * 100 * exposure_time * 0.01
    
    #############################################
    # Convert Wavelength Value to a Pixel Value #
    #      Centered at 1217.67 Angstroms        #
    #############################################
    
    # Conversion factor
    pix_per_ang = col_pixels / ang_range         # Pixels per Angstrom
    
    # Pixel of lyman alpha center wavelength and range of gaussian
    lym_pix = int((lym_alpha_wave - min_x_ang) * pix_per_ang)    # Pixels
    lym_pix_range = int(pix_per_ang)
    
    ###################
    # Create Gaussian #
    ###################
    
    # Create a gaussian the same length as number of columns
    lym_gaussian = create_gaussian(lym_pix_range, lym_pix, lym_std_dev)
    
    ########################
    # Add Padding of Zeros #
    ########################
    
    # Make padding so that gaussian centers at 1217 A
    pad_width_left = lym_pix - (lym_pix_range // 2)
    pad_width_right = col_pixels - pad_width_left - lym_pix_range
    padding = [(pad_width_left, pad_width_right)]
    
    # Create a padded array
    lyman_bg = np.pad(lym_gaussian, pad_width=padding, mode='constant', constant_values=0)
    
    # Multiply by brightness of background
    lyman_bg_array = lyman_bg * lym_photons
    
    # Call convolve function to plot
    lyman_bg, total_counts_per_sec, brightness = (
        convolve_data(lyman_bg_array, slit_data, True, None, None))
    
    return lyman_bg


##############################
# combine background sources #
##############################
def add_background(noisy_2d_spectrum, slit_data):
    """
    Combine background sources onto one spectrum

    Parameters
    ----------
    noisy_2d_spectrum: 2D numpy array
        convolved and normalized SNR / star
        data with added shot noise

    slit_data: 2D numpy array
        slit data file with 760 rows and 1380 cols
        of detector data

    Returns
    -------
    spec_with_background: 2D array
        contains photon data per pixel with
        background noise added

    """
    
    #############################
    # call background functions #
    #############################
    
    dark_bg = gen_dark_bg(noisy_2d_spectrum)
    lyman_bg = gen_lyman_bg(slit_data)
    
    # Add each spectrum together
    spec_with_background = dark_bg + lyman_bg
    
    return spec_with_background
