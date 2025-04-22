###########
# imports #
###########

# Local imports
from calibration_config import noise
global noise

# Third-party Imports
import numpy as np

######################
# noise calculations #
######################
def calculate_noise(photons, scc_photons, waves, exposure_time):
    """
    Calculates signal-to-noise ratio

    Parameters
    ----------
    photons : numpy array
        photon counts corresponding to wave bins

    scc_photons: numpy array
        scc photon counts corresponding to wave bins

    waves : numpy array
        binned wavelengths

    exposure_time: integer
        number of seconds satellite observes star


   Returns
    -------
    snr: numpy array
        signal-to-noise ratio per angstrom

    scc_snr: numpy array
        SCC signal-to-noise ratio per angstrom

    waves: numpy array
        binned wavelengths
    """

    # Empty list to append SNR to
    snr_list = []

    # Empty list ot append 1200 +  wavelengths to
    wave_list = []

    # Empty list to append photon counts to
    photon_count_list = []

    ###################
    # Noise Variables #
    ###################

    # SPRITE noise counts per angstrom
    sprite_noise_per_ang = 0.00006 * exposure_time

    # Noise counts per angstrom
    noise_per_ang = noise * exposure_time

    # Total noise (counts/A)
    total_noise_per_ang = sprite_noise_per_ang + noise_per_ang

    for index in range(len(photons)):

        # Get photon count for each bin
        photon_count = photons[index]

        photon_count_list.append(photon_count)

        # Noise counts
        noise_count = total_noise_per_ang

        # Signal to dark noise ratio
        # Do not take square root of a negative

        if (photon_count + noise_count) <= 0:
            snr = 0

        else:
            snr = photon_count / np.sqrt(photon_count + noise_count)

        # Append SNR to list
        snr_list.append(snr)

        # Append corresponding wavelength to list
        wave_list.append(waves[index])

    snr = np.array(snr_list)
    waves = np.array(wave_list)

    # find the nearest index to 1300
    index_1300 = np.abs(waves - 1300).argmin()

    # signal to noise ratio, average S/NR
    print(f"The S/NR at 1300A = {snr[index_1300]:.2f} for an exposure time of {exposure_time / 1000} kilo-seconds")
    print(f"The average S/NR is: {np.average(snr):.2f} for an exposure time of {exposure_time / 1000} kilo-seconds")

    ###########
    # SCC SNR #
    ###########

    # SUM SCC photon counts for all lambda
    integrated_scc_counts = np.sum(scc_photons)
    # calculate signal to noise for integrated SCC_counts
    wave_nums = len(waves)
    integrated_noise = total_noise_per_ang * wave_nums

    scc_snr = integrated_scc_counts / np.sqrt(integrated_scc_counts + integrated_noise)
    print(f"The total SCC S/NR is: {scc_snr:.2f} for an exposure time of {exposure_time / 1000} kilo-seconds")
    print()

    return snr, scc_snr, waves
