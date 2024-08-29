"""
tools to bin wavelength and flux data and convert to flux to photons
"""

# Third Party Imports
import numpy as np
import astropy.units as u
from astropy.constants import h, c

# Local imports
from sprite_sim_config import *

# Physics constants
h = (h.value * u.joule * u.s).to(u.erg * u.s).value
c = c.value

# Globalize exposure time
global exposure_time


####################
# Binning Function #
####################
def bin_flux_and_wave_data(ang_per_bin, wave_array, flux_array):
    """
    This function takes in wavelength and flux data arrays and
    bins the data

    Parameters
    ---------
    ang_per_bin : integer
        Angstrom range per bin

    wave_array : numpy array
        An array of wavelength values (in angstroms) that
        correspond to specific flux data
        
    flux_array: numpy array
        An array of flux values (ergs/cm^2/s/Angstrom)
        that correspond to specific wavelengths

    Returns
    -------
    avg_bin : numpy array
        1D array with N elements
        array that specifies the bins for
        the wavelength data in angstroms
        
    flux_avg: numpy array
        1D array with N elements
        each element is the average flux per bin
        
    flux_sum : numpy array
        1D array with N elements 
        each element is the total flux sum for each bin
    """

    # Determine min and max wave values
    min_wave = int(np.min(wave_array))
    max_wave = int(np.max(wave_array))

    #####################
    # Initialize Arrays #
    #####################

    # Create array of bin edges 
    bin_edges = np.arange(min_wave, (max_wave + ang_per_bin), ang_per_bin)

    # Determine bin indices that each wave falls into
    bin_indices = np.digitize(wave_array, bin_edges, right=True) - 1

    # Number of bins
    bin_num = len(bin_edges) - 1

    # Ensure bin_indices are within the range of [0, bin_num-1]
    bin_indices = np.clip(bin_indices, 0, bin_num - 1)

    # zeros to store total flux per bin
    flux_sum = np.zeros(bin_num)

    # zeros to add number of flux values per bin
    counts = np.zeros(bin_num)

    #################################
    # Sum Fluxes and Counts per Bin #
    #################################

    np.add.at(flux_sum, bin_indices, flux_array)
    np.add.at(counts, bin_indices, 1)

    ##############################
    # Remove Noise and Negatives #
    ##############################

    # Set negatives to zero
    negatives = flux_sum < 0
    flux_sum = np.where(negatives, 0, flux_sum)

    # Set zeros to next value
    zeros = flux_sum <= 0
    flux_sum = np.where(zeros, np.roll(flux_sum, -1), flux_sum)

    ##################################
    # Calculate Average Flux per Bin #
    ##################################

    flux_avg = np.divide(flux_sum, counts, where=counts != 0)

    # Average bins so that plot is not shifted
    bin_avg = (bin_edges[:-1] + bin_edges[1:]) / 2

    return bin_avg, flux_avg, flux_sum


#####################
# Photon Conversion #
#####################

def convert_flux_to_photons(waves, flux, eff_areas, eff_area_waves):
    """
    converts flux values to photon counts
    ** assumes flux in units of ergs/cm^2/s **

    Parameters
    ----------  
    waves: numpy array
        wavelengths binned per angstrom
        
    flux: 1D numpy array
        fluxes binned according to corresponding wavelengths
        
    eff_areas: numpy array
        SPRITE effective areas for different wavelength values
        
    eff_area_waves: numpy array
        SPRITE wavelength values that correspond to
        different effective areas
    
    Returns
    -------
    photons: numpy array
        photon counts per angstrom corresponding
        to binned wavelength array   
       
    """

    #####################################################
    # Interpolate Data to be in SPRITE Wavelength Range #
    #####################################################

    # Evaluate data at SPRITE wavelengths, make area array same length as binned data
    new_area = np.interp(waves, eff_area_waves, eff_areas)  # cm^2

    #######################
    # Photon Calculations #
    #######################

    # Multiply A-eff * flux to get power
    power = flux * new_area  # Erg/s

    # Determine energy from power and exposure time
    energy = power * exposure_time  # Erg

    # Convert wavelength to units of meters
    wavelength_meter = waves * 1e-10  # Meters

    # Convert from energy to photons
    photon_energy = (h * c) / wavelength_meter  # Energy per photon (ergs)

    # Number of photons
    photons = energy / photon_energy  # Counts

    return photons
