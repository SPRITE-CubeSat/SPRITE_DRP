###########
# imports #
###########

# Third-party Imports
from astropy.io import fits

import numpy as np
import pandas as pd

############################
# aperture/slit constants #
############################

from bin_and_convert_star_data import bin_flux_and_wave_data

from calibration_config import file_path


###############################################
# open fits file and extract wave/photon data #
###############################################

def get_data_from_1_file(file_name, data_base):
    """
    opens and uploads supernova file data

    Parameters
    ----------

    file_name: string
        directory + filename of HST standard star data

    data_base: string
        name of database data was taken from (ie. HST, GALEX, IUE)

    Returns
    -------
    waves: numpy array
        wavelength values (in angstroms) that
        correspond to specific flux data

    flux: numpy array
        An array of flux values (ergs/cm^2/s/Angstrom)
        that correspond to specific wavelengths

    target_name: string
        specifies name of calibration target
    """

    ###############################################
    #                load in data                 #
    ###############################################

    # Open FITs file
    my_fits = fits.open(file_name)

    # Export header info into array
    spec_hdu = my_fits[1]

    if data_base == "HST":

        # Export wavelength and flux data
        wave_data = spec_hdu.data["WAVELENGTH"][0]  # Angstroms
        flux_data = spec_hdu.data["FLUX"][0]  # * (100/6.25)     Erg/cm**2/s/angstrom

        # Get target name (Use HST header)
        with fits.open(file_name) as hdul:

            # data in array format
            header = hdul[0].header
            target_name = header['TARGNAME']

    elif data_base == "IUE":
        wave_data = spec_hdu.data["WAVE"][0]  # Angstroms
        flux_data = spec_hdu.data["FLUX"][0]  # * (100/200)      Erg/cm**2/s/angstrom

        # Get target name (Use HST header)
        with fits.open(file_name) as hdul:

            # data in array format
            # hdul.info()
            header = hdul[0].header
            target_name = header.get('LTARGET', None)  # Returns None if the key is missing
            if not target_name:
                target_name = "Zeta Cas"

    elif data_base == "FUSE":
        wave_data = spec_hdu.data["WAVE"][0]  # Angstroms
        flux_data = spec_hdu.data["FLUX"][0]  # * (100/900)    Erg/cm**2/s/angstrom

        # Get target name (Use HST header)
        with fits.open(file_name) as hdul:

            # data in array format
            # hdul.info()
            header = hdul[0].header
            target_name = header.get('TARGNAME', None)  # Returns None if the key is missing

    elif data_base == "TUES":

        with fits.open(file_name) as hdul:
            wave_data = []  # List to store wavelengths
            flux_data = []  # List to store corresponding flux values
            header = hdul[0].header
            target_name = header["OBJECT"]

            # Loop over each order (HDU 1 to last HDU)
            for i in range(1, len(hdul)):
                data = hdul[i].data  # Access the table data

                # Extract wavelength and flux (modify column name as needed)
                wavelengths = data['WAVELENGTH']  # Wavelength in Angstroms
                fluxes = data['ENERGY_FLUX']  # Photon flux in photons/cm²/s/Å

                wave_data.append(wavelengths)
                flux_data.append(fluxes)

            # Flatten and sort the arrays to get a continuous spectrum
            wave_data = np.concatenate(wave_data)
            flux_data = np.concatenate(flux_data)

            # Sort by wavelength in case orders aren't sequential
            sorted_indices = np.argsort(wave_data)
            wave_data = wave_data[sorted_indices]
            flux_data = flux_data[sorted_indices]

    elif data_base == "HUT":

        # Export wavelength and flux data
        wave_data = spec_hdu.data["WAVE"]  # Angstroms
        flux_data = spec_hdu.data["FLUX"]  # * (100/560) Erg/cm**2/s/angstrom

        # Get target name (Use HST header)
        with fits.open(file_name) as hdul:

            # data in array format
            header = hdul[0].header
            target_name = header['TARGNAME']

    #####################
    # Call Binning Func #
    #####################

    # Bin data per 2 ang (wavelengths)
    waves, flux, sum_flux = bin_flux_and_wave_data(2, wave_data, flux_data)

    ###########################
    # Remove extreme outliers #
    ###########################

    for index in range(len(flux)):
        if flux[index] > (4 * np.mean(flux)):
            flux[index] = flux[index - 1]

    #####################
    # Remove NaN Values #
    #####################

    # Mask is true if value is not NaN
    valid_values = ~np.isnan(flux)
    waves = waves[valid_values]

    if data_base == "IUE":
        flux = sum_flux[valid_values]

    else:
        flux = flux[valid_values]

    return waves, flux, target_name


def get_data_from_2_files(file_fuse, file_iue):

    """
    opens and uploads two data files for a star

    Parameters
    ----------
    file_fuse: string
        directory + filename of FUSE target star data
        for a specific supernova remnant

    file_iue: string
        directory + filename of IUE target star data
        for the same supernova remnant

    Returns
    -------
    waves: numpy array
        wavelength values (in angstroms) that
        correspond to specific flux data

    flux: numpy array
        An array of flux values (ergs/cm^2/s/Angstrom)
        that correspond to specific wavelengths

    target_name: string
        specifies name of calibration target
    """

    ###############################################
    #                load in data                 #
    ###############################################

    #############
    # FUSE FILE #
    #############

    # Open FITs file
    my_fits_fuse = fits.open(file_fuse)

    # Export header info into array
    fuse_spec_hdu = my_fits_fuse[1]

    # Export wavelength and flux data
    fuse_wave_data = fuse_spec_hdu.data["WAVE"][0]  # Angstroms
    fuse_flux_data = fuse_spec_hdu.data["FLUX"][0]  # Erg/cm**2/s/angstrom

    #############
    # IUE FILE #
    #############

    # Open FITs file
    my_fits_iue = fits.open(file_iue)

    # Export header info into array
    iue_spec_hdu = my_fits_iue[1]

    # Export wavelength and flux data
    iue_wave_data = iue_spec_hdu.data["WAVE"][0]  # Angstroms

    iue_flux_data = iue_spec_hdu.data["FLUX"][0]  # Erg/cm**2/s/Angstrom

    ###############
    # Target Name #
    ###############

    # Get target name (use IUE header)
    with fits.open(file_iue) as hdul:
        # data in array format
        # hdul.info()
        header = hdul[0].header
        target_name = header['LTARGET']

    #####################
    # Call Binning Func #
    #####################

    # Bin data per wavelengths
    bins_fuse, avg_flux_fuse, sum_flux_fuse = bin_flux_and_wave_data(2, fuse_wave_data, fuse_flux_data)

    bins_iue, avg_flux_iue, sum_flux_iue = bin_flux_and_wave_data(2, iue_wave_data, iue_flux_data)

    ###########################
    # Remove overlapping bins #
    ###########################

    # Starting IUE wavelength
    min_iue_wave = np.min(bins_iue)

    # Where FUSE wavelengths overlap IUE waves
    overlapping_indices = np.where(bins_fuse >= min_iue_wave)

    # Remove overlapping data
    bins_fuse = np.delete(bins_fuse, overlapping_indices)
    avg_flux_fuse = np.delete(avg_flux_fuse, overlapping_indices)

    ####################
    # Concatenate Data #
    ####################

    # Append together data
    waves = np.concatenate((bins_fuse, bins_iue))
    flux = np.concatenate((avg_flux_fuse, sum_flux_iue))

    #####################
    # Remove NaN Values #
    #####################

    # Mask is true if value is not NaN
    valid_values = ~np.isnan(flux)
    waves = waves[valid_values]
    flux = flux[valid_values]

    return waves, flux, target_name


def get_eff_area_data():
    """
    Opens and uploads effective area data for SPRITE

    Returns
    -------
    sprite_waves: numpy array
        wavelength values (in angstroms) that correspond to specific effective area data for SPRITE

    sprite_scc_waves: numpy array
        wavelength values (in angstroms) that correspond to specific SCC area data for SPRITE

    sprite_eff_area: numpy array
        area values (cm^2) that correspond to SPRITE wavelengths

    sprite_scc_area: numpy array
        SCC area values (cm^2) that correspond to SPRITE wavelengths
    """

    #####################
    # Read in Data File #
    # SPRITE Eff. Area  #
    #####################

    # .dat file with effective area
    eff_area_file = r"/SPRITE_effa.dat"
    scc_eff_area_file = r"/SPRITE_effa_star.dat"

    # Path to file
    eff_area_path = file_path + eff_area_file
    scc_eff_a_path = file_path + scc_eff_area_file

    dat_file = pd.read_csv(eff_area_path, sep="\s+", skiprows=1, header=0,
                           names=["Wavelength", "Grating", "DQE", "M1_LiF", "eLiF", "A_eff", "BEF"])
    scc_dat_file = pd.read_csv(scc_eff_a_path, sep="\s+", skiprows=1, header=0,
                               names=["Wavelength", "Grating", "DQE", "M1_LiF", "eLiF", "A_eff", "BEF", "SCC"])

    sprite_waves = dat_file["Wavelength"]  # Angstrom
    sprite_scc_waves = scc_dat_file["Wavelength"]  # Angstrom
    sprite_eff_area = dat_file["A_eff"]
    sprite_scc_area = scc_dat_file["SCC"]

    return sprite_waves, sprite_scc_waves, sprite_eff_area, sprite_scc_area
