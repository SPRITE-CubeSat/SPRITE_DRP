"""

Tools to extract SNR and slit data

"""

#######################
# Third-party Imports #
#######################

# Pandas
import pandas as pd

# Astropy
from astropy.io import fits

#################
# Local imports #
#################

from bin_and_convert_data import *
from sprite_sim_config import *

# detector constants
global file_path


###############################################
# open fits file and extract wave/photon data #
###############################################
def load_iue_and_fuse_spec(fuse_file, iue_file):
    
    """
    opens and uploads supernova file data

    Parameters
    ----------
    fuse_file: string
        filename of FUSE supernova data
        for a specific supernova remnant

    iue_file: string
        directory + filename of IUE supernova data
        for the same supernova remnant

    Returns
    -------

    waves: numpy array
        An array of wavelength values (in angstroms) that
        correspond to specific flux data for supernova remnant

    flux: numpy array
        An array of flux values (ergs/cm^2/s/Angstrom)
        that correspond to specific wavelengths

    target_name: string
        specifies name of supernova remnant target

    """
    ###############################################
    #                load in data                 #
    ###############################################

    #############
    # FUSE FILE #
    #############
    
    # Open FITs file
    fuse_data = file_path + fuse_file
    my_fits_fuse = fits.open(fuse_data)

    # Export header info into array
    fuse_spec_hdu = my_fits_fuse[1]

    # Export wavelength and flux data        
    fuse_waves = fuse_spec_hdu.data["WAVE"][0]  # Angstroms
    fuse_flux = fuse_spec_hdu.data["FLUX"][0]  # Erg/cm**2/s/angstrom

    #############
    # IUE FILE #
    #############
    
    # Open FITs file
    iue_data = file_path + iue_file
    my_fits_iue = fits.open(iue_data)
    
    # Export header info into array
    iue_spec_hdu = my_fits_iue[1]

    # Export wavelength and flux data        
    iue_waves = iue_spec_hdu.data["WAVE"][0]    # Angstroms

    iue_flux = iue_spec_hdu.data["FLUX"][0]    # Erg/cm**2/s/Angstrom

    ###############
    # Target Name #
    ###############
    
    # Get target name (use FUSE assuming targets are the same)
    with fits.open(fuse_data) as hdul:

        # data in array format
        header = hdul[0].header
        target_name = header['TARGNAME']

    #####################
    # Call Binning Func #
    #####################
    
    # Bin data per wavelengths
    fuse_bins, fuse_avg_flux, fuse_flux_sum = bin_flux_and_wave_data(1, fuse_waves, fuse_flux)

    iue_bins,  iue_avg_flux,  iue_flux_sum = bin_flux_and_wave_data(1, iue_waves, iue_flux)

    # Scale IUE data to better match FUSE SNR data
    iue_flux_sum *= 6
    
    ###########################
    # Remove overlapping bins #
    ###########################
   
    # Starting IUE wavelength
    min_iue_wave = np.min(iue_bins)
    
    # Where FUSE wavelengths overlap IUE waves
    overlapping_indices = np.where(fuse_bins >= min_iue_wave)
    
    # Remove overlapping data
    fuse_bins = np.delete(fuse_bins, overlapping_indices)
    fuse_avg_flux = np.delete(fuse_avg_flux, overlapping_indices)

    ####################
    # Concatenate Data #
    ####################
    
    # Append together data
    waves = np.concatenate((fuse_bins, iue_bins))
    flux = np.concatenate((fuse_avg_flux, iue_flux_sum))
    
    #####################
    # Remove NaN Values #
    #####################
    
    # Mask is true if value is not NaN
    valid_values = ~np.isnan(flux)
    waves = waves[valid_values]
    flux = flux[valid_values]

    return waves, flux, target_name


##############################################
# open effective area file and extract areas #
##############################################
def load_eff_area(eff_area_file):
    
    """
    loads in effective area data for corresponding wavelengths
    (area in units of cm^2)

    Parameters
    ----------
    eff_area_file: string
        filename of SPRITE effective area data file

    Returns
    -------
    eff_areas: numpy array
        SPRITE effective areas for different wavelength values

    eff_area_waves: numpy array
        SPRITE wavelength values that correspond to
        different effective areas
   """

    ####################
    # Read in Data ile #
    # SPRITE Eff. Area  #
    #####################
    
    eff_area = file_path + eff_area_file
    dat_file = pd.read_csv(
        eff_area, sep=r"\s+", skiprows=1, header=0,
        names=["Wavelength", "Grating", "DQE", "M1_LiF", "eLiF", "A_eff", "BEF"])
    
    eff_area_waves = dat_file["Wavelength"]
    eff_areas = dat_file["A_eff"]
    
    return eff_areas, eff_area_waves


#######################
# Open Slit Data File #
#######################
def open_slit_file(slit_file):
    
    """
    Opens binary slit data file and stores data as a 2D array

    Parameters
    ----------
    slit_file: string
        contains directory and filename of binary
        slit data file

    Returns
    -------
    good_data: 2D numpy array
        data values that store binary detector info
        with all other values removed
    """
    ##################
    # FITS File Data #
    ##################

    slit_file = file_path + slit_file
    with fits.open(slit_file) as hdu_list:
        
        # Data in array format
        data = hdu_list[0].data
        
    # Only look at actual data
    good_data = data[855:1615, 320:1700]
    
    return good_data


########################
# Load brightness data #
########################
def load_brightness_profile(profile):
    
    """
    Parameters:
    -----------
    profile: string
        Name of brightness profile file

    Returns:
    --------
    e_per_s: numpy array
        electron per second count for a given coordinate
        
    original_x_coords: numpy array
        positions corresponding to electron per second counts
    """

    # Path to brightness profile
    brightness_path = file_path + profile
    
    # Get brightness profile from data
    data_file = pd.read_csv(brightness_path, sep=r"\s+", header=None, names=["location", "e/s"])
    
    # Electrons/second
    e_per_s = data_file["e/s"]
    
    # X-coordinates
    original_x_coords = np.array(data_file["location"])
    
    return e_per_s, original_x_coords
