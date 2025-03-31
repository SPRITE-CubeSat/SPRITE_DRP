from get_star_data import *
from bin_and_convert_star_data import *
from noise_and_calculations import *
from visualize_data import *
global file_path


######################
# call all functions #
######################
def call_functions_one_database(file_name, exposure_times, data_base):
    """
    This function is called if only one datafile is used
    (ex: HST or IUE data)

    Parameters
    ---------
    file_name : string
        name of the flux data file being analyzed

    exposure_times: list
        exposure time in seconds

    data_base: string
        mission where datafile was obtained

    Returns
    -------
    None
    """

    # Empty lists to append data for each exposure times
    pixel_list = []
    flux_list = []
    photon_per_ang_list = []
    photon_per_pix_list = []
    snr_list = []
    wave_list = []

    for time in exposure_times:
        # Get star data
        waves, flux, target_name = get_data_from_1_file(file_name, data_base)

        # Get effective area data
        sprite_waves, sprite_scc_waves, sprite_eff_area, sprite_scc_area = get_eff_area_data()

        # Convert flux/A to photons/A
        photons, scc_photons = flux_to_photons(sprite_waves, sprite_scc_waves, waves, sprite_eff_area, sprite_scc_area,
                                               flux, time)

        # Calculate SNR values and photon counts
        snr, scc_snr, waves = calculate_noise(photons, scc_photons, waves, time)

        # Get pixel values
        pixels, photons_per_pix, photons_per_ang, flux_per_ang, clipped_waves, clipped_snr = clip_waves_and_convert_to_pix(
            waves, flux, photons, snr)

        # append to lists
        flux_list.append(flux_per_ang)
        pixel_list.append(pixels)
        photon_per_ang_list.append(photons_per_ang)
        photon_per_pix_list.append(photons_per_pix)
        snr_list.append(clipped_snr)
        wave_list.append(clipped_waves)

    # average flux
    avg_flux = np.mean(flux_list[0])
    min_wave = np.min(wave_list[0])
    max_wave = np.max(wave_list[0])

    print(f"The average flux from {min_wave} - {max_wave} Angstroms is {avg_flux:.1g}")

    #########
    # Plots #
    #########
    plot_photons_vs_pix(pixel_list, photon_per_pix_list, exposure_times, target_name)

    plot_snr(snr_list, wave_list, exposure_times, target_name)

    plot_photons_vs_waves(wave_list, photon_per_ang_list, exposure_times, target_name)

    plot_flux_vs_waves(wave_list[0], flux_list[0], target_name)


def call_functions_two_database(file_fuse, file_iue, exposure_times):
    """
    This function is used instead of call_functions_one_database if using more
    than one file (ex: FUSE and IUE data)

    Parameters
    ----------

    file_fuse: string
        directory + filename of FUSE supernova data
        for a specific supernova remnant

    file_iue: string
        directory + filename of IUE supernova data
        for the same supernova remnant

    exposure_times: list
        specifies exposure times of satellite in seconds

    Returns
    -------
    None
    """

    # Empty lists to append data for dif exposure times

    flux_list = []
    pixel_list = []
    photon_per_ang_list = []
    photon_per_pix_list = []
    snr_list = []
    wave_list = []

    for time in exposure_times:
        # Get star data
        waves, flux, target_name = get_data_from_2_files(file_fuse, file_iue)

        # Get eff_area data
        sprite_waves, sprite_scc_waves, sprite_eff_area, sprite_scc_area = get_eff_area_data()

        # Convert flux/A to photons/A
        photons, scc_photons = flux_to_photons(sprite_waves, sprite_scc_waves, waves, sprite_eff_area, sprite_scc_area, flux, time)

        # Calculate SNR values and photon counts
        snr, scc_snr, waves = calculate_noise(photons, scc_photons, waves, time)

        # Get pixel values
        pixels, photons_per_pix, photons_per_ang, flux_per_ang, clipped_waves, clipped_snr = clip_waves_and_convert_to_pix(
            waves, flux, photons, snr)

        # append to lists
        flux_list.append(flux_per_ang)
        pixel_list.append(pixels)
        photon_per_ang_list.append(photons_per_ang)
        photon_per_pix_list.append(photons_per_pix)
        snr_list.append(clipped_snr)
        wave_list.append(clipped_waves)

    # average flux
    avg_flux = np.mean(flux_list[0])
    min_wave = np.min(wave_list[0])
    max_wave = np.max(wave_list[0])

    print(f"The average flux from {min_wave} - {max_wave} Angstroms is {avg_flux:.1g}")
    #########
    # Plots #
    #########

    # plot_photons_vs_pix(pixel_list, photon_list, exposure_times, target_name)

    plot_snr(snr_list, wave_list, exposure_times, target_name)

    plot_photons_vs_pix(pixel_list, photon_per_pix_list, exposure_times, target_name)

    plot_photons_vs_waves(wave_list, photon_per_ang_list, exposure_times, target_name)

    plot_flux_vs_waves(wave_list[0], flux_list[0], target_name)
