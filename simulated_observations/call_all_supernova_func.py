"""
Tool to combine all functions for end result of 2D and 1D spectral plots
"""

###########
# imports #
###########

# Local Imports
from load_supernova_data import *
from create_2d_supernova_spec import *
from supernova_spec_noise import *
from plot_supernova_spectra import *
from reconstruct_supernova_image import *

global file_path, fig_path, exposure_time


######################
# Call all functions #
######################
def get_all_supernova_data(slit_file, fuse_file, iue_file, eff_area_file, profile_name, number):
   
    """
    Calls all functions except plotting function
    gets data for all three targets

    Parameters
    ----------       
    slit_file: string
        specifies directory and filename of
        binary slit data file
        
    fuse_file: string
        directory and filename of FUSE data
        for a specific target
        
    iue_file: string
        directory and filename of IUE
        data for a specific target
    
    eff_area_file: string
        directory and filename of SPRITE
        effective area data file
        
    profile_name: string
        name of SNR to extract profile from
        
    number: integer
        number of profiles to loop through

    Returns
    -------
    None: purpose is to call all of the above functions

    """
    
    # Variable to append 2D brightness arrays
    brightness_arrays = []
    
    # Loads in slit data as a 2D array
    slit_data = open_slit_file(slit_file)

    # Loads SNR data into 1D wave and flux arrays
    waves, flux, target_name = load_iue_and_fuse_spec(fuse_file, iue_file)
    
    # Loads SPRITE effective area data into 1D arrays of area and corresponding wavelength
    area, area_waves = load_eff_area(eff_area_file)
    
    # Converts flux array to photon array
    photons = convert_flux_to_photons(waves, flux, area, area_waves)

    # Converts the SN photon counts per wavelength to photon counts per pixel 
    # (1380 pixels in length)
    photon_per_pix = interp_waves_to_pixels(waves, photons)

# If no brightness profile (or a calibration star)
    if profile_name is None:

        # Convolve 1D SNR photon data onto 2D SPRITE slit
        spec_2d, total_counts_per_sec, brightness = convolve_data(photon_per_pix, slit_data, False, None, None)

        # Generate 2D spectrum with shot noise and background
        spec_2d_with_background = generate_noisy_image(spec_2d, slit_data)

        # Plot 1D simulated spec
        plot_1d_spectrum(waves, photons, total_counts_per_sec, target_name)

        # Plot noisy 2D spectrum
        plot_noisy_2d_spec(spec_2d_with_background, number, target_name)

# If SNR with brightness profiles
    else:
        for num in range(number):

            profile = profile_name + f"_{num + 1}.dat"

            # Open a specific brightness profile, get data into arrays
            electron_per_sec, x_coord = load_brightness_profile(profile)

            # Convolve 1D SNR photon data onto 2D SPRITE slit
            spec_2d, total_counts_per_sec, brightness = (
                convolve_data(photon_per_pix, slit_data, False, electron_per_sec, x_coord))

            # Plot the 1D spectrum for the first profile only
            if num == 0:
                plot_1d_spectrum(waves, photons, total_counts_per_sec, target_name)

            # Generate 2D spectrum with shot noise and background
            spec_2d_with_background = generate_noisy_image(spec_2d, slit_data)

            # Plot noisy 2D spectrum
            plot_noisy_2d_spec(spec_2d_with_background, num, target_name)

            # Make 1D brightness profiles 2D (40 x 760 pixels)
            interpolated_counts_2d = convert_brightness_profile_to_2d_array(electron_per_sec)

            # Append profiles to a list
            brightness_arrays.append(interpolated_counts_2d)

        # Stack brightness profiles together and create image
        flipped_array = stack_slit_data_to_reconstruct_supernova_image(brightness_arrays)

        # Plot reconstructed image
        plot_reconstructed_images(flipped_array, target_name)
