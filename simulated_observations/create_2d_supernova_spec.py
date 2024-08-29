"""

Tools to create a 2D array of SN remnant data convolved onto SPRITE slit

"""

#######################
# Third-party Imports #
#######################

# Numpy
import numpy as np

# Scipy
from scipy.signal import convolve, gaussian
from scipy.interpolate import interp1d

#################
# Local imports #
#################

from sprite_sim_config import *

# detector constants
global convolve_constants, exposure_time, quantum_efficiency, supernova_height


#################################################
# Convert Photons Per Wave to Photons Per Pixel #
#################################################

def interp_waves_to_pixels(waves, photons):
    """
    Converts the supernova photon counts per wavelength
    to photon counts per pixel (1380 pixels in length)

    Parameters
    ----------   
    waves: 1D numpy array
        Wavelengths for SNR / star data (binned per A)
        
    photons: 1D numpy array
        Photons corresponding to each wavelength value


    Returns
    -------
    photon_per_pix: 1D numpy array
        Contains photons counts corresponding to each column pixel
    """

    ##########################
    # Min/Max Angstrom Range #
    ##########################

    # Pixels per angstrom
    pix_per_ang = col_pixels / ang_range

    # Min and max wavelength of our SN data
    min_wave = np.min(waves)
    max_wave = np.max(waves)

    ######################################################
    #    check that waves are within aperture bounds    #
    # if not, reset min and max to aperture min and max #
    ######################################################

    # If minimum wave is less than minimum aperture wave
    if min_wave < min_x_ang:
        # Difference in wavelengths
        index1 = int(min_x_ang - min_wave)

        # New min wave = aperture min wave
        min_wave = min_x_ang

        # Get rid of data that will not be displayed
        photons = photons[index1:]
        waves = waves[index1:]

    if max_wave > max_x_ang:
        # Difference in wavelengths
        index2 = int(max_wave - max_x_ang)

        # New min wave = aperture min wave
        max_wave = max_x_ang

        # Get rid of data that will not be displayed
        photons = photons[:(-index2)]
        waves = waves[:(-index2)]

    ################################
    # Angstrom -> Pixel Conversion #
    ################################

    # Min and max wavelength converted to pixel
    min_pix = int((min_wave - min_x_ang) * pix_per_ang)
    max_pix = int((max_wave - min_x_ang) * pix_per_ang)

    # SN pixel range of data
    supernova_pix_range = max_pix - min_pix

    # Create a wavelength array that is the same length as the number of pixels our data ranges
    wave_array_pix_len = np.linspace(min_wave, max_wave, supernova_pix_range)

    # Interpolate photon counts to correspond to each pixel in 2048-long pixel array
    supernova_photons = np.interp(wave_array_pix_len, waves, photons)

    # Add new data to array of zeros
    photon_per_pix = [0] * col_pixels
    photon_per_pix[min_pix:max_pix] = supernova_photons

    return photon_per_pix


################################
# Create a Normalized Gaussian #
################################

def create_gaussian(array_length, center, std_dev):
    """
    Creates a gaussian distribution for a specified array length, mean, and std_dev
    
    Parameters:
    ----------
    array_length: integer
        number of array values that will make up gaussian
        
    center: integer
        mean / center of gaussian
    
    std_dev: integer
        standard deviation of gaussian
    
    Returns
    -------
    gaussian: numpy array
        normalized gaussian distribution created from the above parameters
    
    """

    # Generate an array of x values
    x = np.linspace(0, array_length - 1, array_length)

    # Create a gaussian 
    gaussian = np.exp(-0.5 * ((x - center) / std_dev) ** 2)

    # Normalize the Gaussian array
    gaussian /= gaussian.sum()

    return gaussian


####################################################
# Make a 1D Array of SNR Brightness Along the Slit #
####################################################

def extract_brightness_data(e_per_s, original_x_coords):
    """
    Extracts brightness data along the 1800" slit and converts
    to photons per pixel for 760 row pixels
    
    Parameters:
    -----------
    e_per_s: numpy array
        electron per second count for a given coordinate
        
    original_x_coords: numpy array
        positions corresponding to electron per second counts

    Returns:
    --------
    new_brightness_prof: numpy array
        brightness profile of SNR that is same 
        length as number of row pixels
        
    new_x_coords: numpy array
        row numbers that supernova remnant spans
        
    padded_array: numpy array
        includes SNR data and padding of zeros
        to fill all rows (760 pixels) on image 
        
    """

    # Convert from electrons/s to photons (assume quantum efficiency of 0.5)
    photons = (e_per_s * exposure_time) / quantum_efficiency

    # Remove noise
    outliers = np.where(photons > (4 * np.median(photons)))[0]
    photons[outliers] = np.mean(photons)

    # Determine relative size of SNR to slit height of 1800"
    ratio = supernova_height / 1800

    # Determine how many pixels SNR would span in slit
    supernova_row_pix = int(ratio * row_pixels)

    # New x coordinates to match SNR pixels length
    x_length = len(photons)
    new_x_coords = np.linspace(1, x_length, supernova_row_pix)

    # Interpolate to make brightness profile same length as SNR pixel height
    interp_func = interp1d(original_x_coords, photons, kind='linear')

    # Resampled brightness profile
    new_brightness_prof = interp_func(new_x_coords)

    # Add padding
    # Make padding uneven to simulate real SNR view
    pad_width_left = (row_pixels - supernova_row_pix) // 3
    pad_width_right = row_pixels - supernova_row_pix - pad_width_left
    padding = [(pad_width_left, pad_width_right)]

    # Create a padded array
    padded_array = np.pad(new_brightness_prof, pad_width=padding, mode='constant', constant_values=0)

    return new_brightness_prof, new_x_coords, padded_array


##################################
# Convolve Photon Data onto Slit #
##################################

def convolve_data(photon_per_pix, slit_data, noise=False, e_per_s=None, original_x_coords=None):
    """
    Normalizes SNR / star photon data and convolve it onto the SPRITE slit

    Parameters
    ----------  
    photon_per_pix: 1D numpy array
        Column pixels containing supernova photon data
        
    slit_data: 2D numpy array
        slit data file with 760 rows and 1380 cols
        of detector data

    noise: boolean
        determines whether geo-coronal lyman alpha noise
        is being convolved onto the slit or 
        regular data

    e_per_s: numpy array OR None value
        electron per second count for a given coordinate
        OR
        None if lyA data
        
    original_x_coords: numpy array
        positions corresponding to electron per second counts
        (in brightness profile)
        OR
        None if lyA data

    Returns
    -------
    spec2D: convolved 2D spectrum
        
    max_count_per_sec_2D: float
        maximum photon count per pixel for 2D spectrum
        divided by total exposure time
    
    total_counts_per_sec: float
        total number of photon counts divided by total 
        exposure time for the 1D spectrum
        
    brightness: float
        describes increase/decrease in star brightness
        randomly determined
        
    """

    ####################################
    # Determine Slit Loc. and Binarize #
    ####################################

    # Slit location
    slit_image = slit_data[:, 1040:1080]

    # Convert to boolean mask
    slit_med = np.median(slit_image)
    slit_mask = slit_image > (3 * slit_med)

    # Convert to binary mask
    slit_mask = slit_mask.astype(np.float64)

    # Determine total counts/exposure time
    total_counts = np.sum(photon_per_pix)
    total_counts_per_sec = total_counts / exposure_time

    ###############################
    # GEO-CORONAL LYM ALPHA NOISE #
    ###############################

    # Stars have high photon count - check if data corresponds to that of a star
    if noise:

        # Lyman Alpha Noise = do nothing
        brightness = 1

        for col in range(40):
            slit_mask[:, col] *= 1

    #############
    # STAR DATA #
    #############
    elif total_counts_per_sec > 400:

        # Do not change star brightness
        brightness = 1
        photon_per_pix = np.array(photon_per_pix) * brightness

        # Create a gaussian around star with point-source function 
        # Assume 14 arcsecond point source
        star_gaussian = gaussian(row_pixels, pix_std_dev)

        # Apply shifted star gaussian to each column of data in slit
        for col in range(40):
            slit_mask[:, col] *= star_gaussian

    ############
    # SNR DATA #
    ############

    else:
        # No change in brightness
        brightness = 1

        # Brightness_curve  = extract_brightness_profile (exposure_time)
        bp, x, brightness_curve = extract_brightness_data(e_per_s, original_x_coords)

        # Apply SN light curve to each column of data
        for col in range(40):
            slit_mask[:, col] *= brightness_curve

    #######################
    # Normalize Slit Data #
    #######################

    # Add up total photon intensity of slit
    slit_sum = np.sum(slit_mask)

    # Normalize total intensity to equal 1
    normalized_slit = slit_mask / slit_sum

    ############
    # Convolve #
    ############

    # Make photon data a 2D array    
    photon_per_pix = np.array([photon_per_pix])

    # Convolve data
    spec2d = convolve(normalized_slit, photon_per_pix, mode="full")

    # Remove extra 40 pixels of data
    spec2d = spec2d[:, 20:-19]

    # Remove negatives
    spec2d[spec2d < 0] = 0

    if not noise:
        # Check if sum of photons is the same for 2D and 1D spectrum
        print(f"The total photon count from the 1D spectrum is: {np.sum(photon_per_pix):.1f}")
        print(f"The total photon count from the 2D spectrum is: {np.sum(spec2d):.1f}")

    return spec2d, total_counts_per_sec, brightness
