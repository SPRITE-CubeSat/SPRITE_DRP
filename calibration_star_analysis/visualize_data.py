import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

from calibration_config import *

global min_x_ang, max_x_ang, figure_path


#########################################
# Plot Signal-Noise Ratio VS Wavelength #
#########################################
def plot_snr(snr, waves, exposure_times, target_name):
    """
    Plot signal-to-noise ratio at each wavelength

    Parameters
    ----------
    snr: nested list
        signal-to-noise ratio per binned wave array for
        each exposure time

    waves: nested list
        binned wavelengths for each exposure time (Angstroms)

    exposure_times: list
        exposure times in seconds

    target_name: string
        name of target that data is extracted from

    Returns
    -------
    None

    """

    pixels = np.arange(0, 2040, 1)
    goal_min_snr = [10] * 2040

    #############
    # S/NR plot #
    #############

    # Initialize figure
    plt.figure(figsize=(9, 6))
    plt.rcParams["font.family"] = "sans serif"
    plt.rcParams["font.size"] = 18

    color = ["#9aaded", "#440154"]

    for index in range(len(exposure_times)):
        # Plot SNR vs Wavelength for each exposure time
        plt.plot(waves[index], snr[index], color[index], lw=2,
                 label=f"Exposure Time = {exposure_times[index]:.0f} seconds")

        # Labels and legend
        plt.xlabel(r"Wavelength $(\AA)$", fontsize=20)
        plt.ylabel("S/N", fontsize=20)
        plt.title(f"Signal to Noise Ratio for {target_name}", fontsize=24)
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

        # Axes limits
        plt.xlim(min_x_ang, max_x_ang)

        # Get the current axes
        ax = plt.gca()

        # Set y-axis to scientific notation
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))

        # Set scientific notation limits
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 3))

    # Plot minimum SNR goal
    plt.plot(pixels, goal_min_snr, color="mediumorchid", linestyle="--", lw=3, label="S/N = 10")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

    # Save figure to file
    plt.savefig(figure_path + f"snr_per_wave_{target_name}.png", bbox_inches="tight")
    # Show figure
    plt.show()


##########################
# Plot Photons per Pixel #
##########################
def plot_photons_vs_pix(new_pixels, padded_photons, exposure_times, target_name):
    """
    Plot photon counts at each pixel

    Parameters
    ----------
    new_pixels: numpy array
       pixel values

    padded_photons: numpy array
        photon counts at different pixels

    exposure_times: list
        exposure times in seconds

    target_name: string
        name of target that data is extracted from

    Returns
    -------
    N/A

    (purpose of function is to plot data)
    """

    ######################
    # Plot Flux vs Waves #
    ######################

    # Initialize figure
    plt.figure(figsize=(9, 6))
    plt.rcParams["font.family"] = "sans serif"
    plt.rcParams["font.size"] = 18

    color = ["#33c46d", "#083f05"]

    for index in range(len(exposure_times)):
        # Plot of wavelength vs flux for FITs file
        plt.plot(new_pixels[index], padded_photons[index], color[index], lw=2,
                 label=f"Exposure time = {exposure_times[index]:.0f} seconds")

        # Axes labels and title
        plt.xlabel("Detector Pixels", fontsize=20)
        plt.ylabel("Photon Counts", fontsize=20)
        plt.title(f"Photon Counts per Pixel for Target: {target_name}", fontsize=24)
        plt.xlim(0, 2040)
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

        # Get the current axes
        ax = plt.gca()

        # Set the formatter for the y-axis to scientific notation
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))

        # Optionally, set the limits for scientific notation
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 3))

    # Save figure to file
    plt.savefig(figure_path + f"photon_per_pix_{target_name}.png", bbox_inches="tight")

    plt.show()


##############################
# Plot Photons VS Wavelength #
##############################
def plot_photons_vs_waves(waves, photons, exposure_times, target_name):
    """
    Plot photon counts at each wavelength

    Parameters
    ----------
    waves: numpy array
        binned wavelengths for each exposure time (Angstroms)

    photons: numpy array
        photon counts per wave value

    exposure_times: list
        exposure times in seconds

    target_name: string
        name of target that data is extracted from

    Returns
    -------
    N/A

    (purpose of function is to plot data)
    """

    ######################
    # Plot Flux vs Waves #
    ######################

    # Initialize figure
    plt.figure(figsize=(9, 6))
    plt.rcParams["font.family"] = "sans serif"
    plt.rcParams["font.size"] = 18

    color = ["#52c4c4", "#2d5abc"]

    for index in range(len(exposure_times)):
        # Plot of wavelength vs flux for FITs file
        plt.plot(waves[index], photons[index], color[index], lw=2,
                 label=f"Exposure time = {exposure_times[index]:.0f} seconds")

        # Axes labels and title
        plt.xlabel(r"Wavelength $(\AA)$", fontsize=20)
        plt.ylabel("Photon Counts", fontsize=20)
        plt.title(f"Photon Counts per Angstrom for Target: {target_name}", fontsize=24)
        plt.xlim(min_x_ang, max_x_ang)
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

        # Get the current axes
        ax = plt.gca()

        # Set the formatter for the y-axis to scientific notation
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))

        # Optionally, set the limits for scientific notation
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 3))

    # Save figure to file
    plt.savefig(figure_path + f"photon_per_wave_{target_name}.png", bbox_inches="tight")
    plt.show()


###########################
# Plot Flux VS Wavelength #
###########################
def plot_flux_vs_waves(waves, flux, target_name):
    """
    Plot flux at each wavelength
    Plot min/max flux values for comparison

    Parameters
    ----------
    waves: numpy array
        binned wavelengths for each exposure time (Angstroms)

    flux: numpy array
        flux at different wavelengths

    target_name: string
        name of target that data is extracted from
    """

    ######################
    # Plot Flux vs Waves #
    ######################

    # initialize figure
    plt.figure(figsize=(9, 6))
    plt.rcParams["font.family"] = "sans serif"
    plt.rcParams["font.size"] = 18

    # plot waves vs flux
    plt.plot(waves, flux, color="red", lw=2)

    # Axes labels and title
    plt.title(f"Flux per Angstrom for Target: {target_name}", fontsize=24)
    plt.xlabel(r"Wavelength $(\AA)$", fontsize=20)
    plt.ylabel("Flux", fontsize=20)
    plt.xlim(min_x_ang, max_x_ang)
    plt.yscale("log")

    plt.axhline(5e-11, color="red", linestyle="--", label=r"flux/$\AA$ = 5e-11")
    plt.axhline(1e-10, color="orange", linestyle="--", label=r"flux/$\AA$ = 1e-10")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

    # Save figure to file
    plt.savefig(figure_path + f"flux_per_wave_{target_name}.png", bbox_inches="tight")
    plt.show()
