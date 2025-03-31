###########
# imports #
###########

# Third-party Imports
import numpy as np
import astropy.units as u
from astropy.constants import h, c

# Physics constants
h = (h.value * u.joule * u.s).to(u.erg * u.s).value
c = c.value

# Local Imports
from calibration_config import *
global row_pixels, col_pixels, min_x_ang, max_x_ang, ang_range


####################
# Binning Function #
####################
def bin_flux_and_wave_data(ang_per_bin, wave_array, flux_array):
    """
    Bins wavelength and flux data for a given angstrom range

    Parameters
    ---------
    ang_per_bin : integer
        Angstrom range per bin

    wave_array : numpy array
        wavelength values (Angstrom)

    flux_array: numpy array
        flux values (ergs/cm^2/s/Angstrom) that correspond to specific wavelengths

    Returns
    -------
    bin_avg : numpy array
        average wavelength per bin

    flux_avg: numpy array
        average flux per bin (erg/cm^2/s/Angstrom)

    flux_sum : numpy array
        the total flux for each bin
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
    N = len(bin_edges) - 1

    # Ensure bin_indices are within the range of [0, N-1]
    bin_indices = np.clip(bin_indices, 0, N - 1)

    # N zeros to store total flux per bin
    flux_sum = np.zeros(N)

    # N zeros to add number of flux values per bin
    counts = np.zeros(N)

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
def flux_to_photons(sprite_waves, sprite_scc_waves, waves, sprite_eff_area, sprite_scc_area, flux, exposure_time):
    """
    converts flux values to photon counts
    ** assumes flux in units of ergs/cm^2/s **

    Parameters
    ----------
    sprite_waves: numpy array
        SPRITE wavelengths corresponding to effective area data

    sprite_scc_waves: numpy array
        SPRITE wavelengths corresponding to SCC effective area data

    waves: numpy array
        binned wavelength values

    sprite_eff_area: numpy array
        SPRITE effective area value per wavelength (cm^2)

    sprite_scc_area
        SPRITE SCC effective area per wavelength (cm^2)

    flux: numpy array
        fluxes binned according to corresponding wave values

    exposure_time: integer
        number of seconds satellite observes star

    Returns
    -------
    photons: numpy array
        photon counts per wavelength for detector

    scc_photons: numpy array
        photon counts per wavelength for SCC detector
    """

    # Evaluate data at SPRITE wavelengths, make area array same length as binned data
    area_new = np.interp(waves, sprite_waves, sprite_eff_area)  # cm^2
    scc_area_new = np.interp(waves, sprite_scc_waves, sprite_scc_area)  # cm^2

    #######################
    #  Photon Calculations #
    #######################

    # Multiply A-eff * flux to get power
    power = flux * area_new  # erg/s
    scc_power = flux * scc_area_new

    # Determine energy from power and exposure time
    energy = power * exposure_time  # erg
    scc_energy = scc_power * exposure_time

    # convert wavelength to units of meters
    wavelength_meter = waves * 1e-10

    # Convert from energy to energy per photon
    photon_energy = (h * c) / wavelength_meter  # erg

    # Number of photons
    photons = energy / photon_energy
    scc_photons = scc_energy / photon_energy

    return photons, scc_photons


def clip_waves_and_convert_to_pix(waves, flux, photons, snr):
    """
    Converts wavelength values to pixels on detector


    Parameters
    ----------
    waves: numpy array
        binned wavelength values

    flux: numpy array
        fluxes binned according to corresponding wave values

    photons: numpy array
        photon counts per wavelength

    snr: numpy array
        signal-to-noise ratio per binned wavelength

    Returns
    -------
    pixels: numpy array
        pixels across entire detector

    photons_per_pix: numpy array
        number of photon counts corresponding to each pixel value

    photons: numpy array
        photon counts per wave bin clipped for SPRITE wavelengths

    flux: numpy array
        flux per wave bin clipped for SPRITE wavelengths

    waves: numpy array
        waves per bin clipped for SPRITE wavelengths

    snr: numpy array
        signal-to-noise ratio per wave bin clipped for SPRITE wavelengths

    """

    ##########################
    # Min/Max Angstrom Range #
    ##########################

    # Pixels per angstrom
    # pix_per_ang = col_pixels / ang_range

    # Array of pixel values
    pixels = np.arange(2040)

    # Min and max wavelength of our SN data
    star_min_wave = np.min(waves)
    star_max_wave = np.max(waves)

    ######################################################
    #    Check that waves are within aperture bounds    #
    # If not, reset min and max to aperture min and max #
    ######################################################

    # If minimum wave is less than minimum aperture wave
    if star_min_wave < min_x_ang:
        # Difference in wavelengths
        index1 = int((min_x_ang - star_min_wave) / 2)

        # New min wave = aperture min wave
        star_min_wave = min_x_ang

        # Get rid of data that will not be displayed
        snr = snr[index1:]
        photons = photons[index1:]
        flux = flux[index1:]
        waves = waves[index1:]

    if star_max_wave > max_x_ang:
        # Difference in wavelengths
        index2 = int(np.ceil((star_max_wave - max_x_ang) / 1))

        # New max wave = aperture max wave
        star_max_wave = max_x_ang

        # Get rid of data that will not be displayed
        snr = snr[:(-(index2 + 1))]
        photons = photons[:(-(index2 + 1))]
        flux = flux[:(-(index2 + 1))]
        waves = waves[:(-(index2 + 1))]

    ################################
    # Angstrom -> Pixel Conversion #
    ################################

    # Min and max wavelength converted to pixel
    # Minimum pixel must be at least 320
    star_min_pix = int((star_min_wave - min_x_ang)) + 320
    star_max_pix = int((star_max_wave - min_x_ang)) + 320

    # Star pixel range of data
    star_pix_range = star_max_pix - star_min_pix

    # Create an array that is the same length as the number of pixels our data ranges
    wave_array_pix_len = np.linspace(star_min_wave, star_max_wave, star_pix_range)

    # Interpolate photon counts to correspond to each pixel in 2048-long pixel array
    star_photons = np.interp(wave_array_pix_len, waves, photons)

    # Add new data to array of zeros
    photon_per_pix = np.zeros(2040)
    photon_per_pix[star_min_pix:star_max_pix] = star_photons

    return pixels, photon_per_pix, photons, flux, waves, snr
