"""
Tools to plot 1D and 2D SNR spectrum and reconstructed SNR image
"""

#######################
# Third-party Imports #
#######################

# Matplotlib / Plotting
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Numpy
import numpy as np

#################
# Local Imports #
#################

from sprite_sim_config import *

global max_x_ang, min_x_ang, fig_path, exposure_time


#####################
# Plot 1D spectrum  #
#####################
def plot_1d_spectrum(waves, photons, total_counts_per_sec, target_name):
    """
    Plot wavelength vs photon count
    (1D plot)

    Parameters
    ----------
    waves: numpy array
       average wavelength value corresponding to photon count array

    photons: numpy array
        photon counts at different wavelengths

    target_name: string
        specifies name of SN remnant target

    total_counts_per_sec: int
        total photon counts per second

    Returns
    -------
    N/A

    (purpose of function is to plot data on a 1D spectrum)
    """

    ##################
    # Plot Constants #
    ##################

    # Determine bounds from wavelengths
    lower_bound = min_x_ang - 3
    upper_bound = max_x_ang + 3

    # Initialize figure
    fig, ax = plt.subplots(figsize=(14, 9))

    #########
    # Plots #
    #########

    # Plot of wavelength vs photon count for supernova file
    plt.plot(waves, photons, color="#0b9ba8")

    # Axes labels
    ax.set_xlabel(r"Wavelength ($\mathring{A}$)", fontsize=30)
    ax.set_ylabel("Photons (counts)", fontsize=30)
    ax.set_title(f"Simulated 1D Spectrum for {target_name}", fontsize=30)
    ax.tick_params(axis='both', which='major', labelsize=24)
    ax.set_xlim(lower_bound, upper_bound)

    # Determine if data is for a star
    if total_counts_per_sec >= 400:
        ax.set_ylim(0, 1100)
    else:
        ax.set_ylim(0, 400)

    # Save figure to file
    plt.savefig(fig_path + f"Wavelength_VS_Photons_{target_name}.png", bbox_inches="tight")

    # Show plot
    plt.show()


##############################
#      plot 2D spectrum      #
#  with background and noise #
##############################

def plot_noisy_2d_spec(spec_2d_with_background, num, target_name):
    """
    Parameters
    ----------
    spec_2d_with_background: 2D numpy array
        contains SNR and star photon data with
        noise and background

    num: integer
        specifies what slit profile we are looking at

    target_name: string
        specifies name of SN target

    Returns
    -------
    None: purpose is to plot noisy 2D spectrum
    """

    # profile number
    num += 1

    # Initialize figure
    fig, ax = plt.subplots(figsize=(14, 9))

    # Display the image using a magma colormap
    image = plt.imshow(spec_2d_with_background,
                       extent=[min_x_ang, max_x_ang, 23, -23], cmap="viridis", vmin=0, vmax=20, aspect="auto")

    # Title and axes
    ax.set_title(f"Simulated 2D Spectrum for {target_name} "
                 f"(t = {(exposure_time / 1000):.0f}ks)\nProfile #{num}", fontsize=30)
    ax.set_xlabel(r'Wavelength ($\mathring{A}$)', fontsize=30)
    ax.set_ylabel('Arc-minutes', fontsize=30)

    # Increase number of x-tick marks to show aperture is lined up at 1600 A
    ax.set_xticks(np.arange(1000, 1801, 100))
    ax.tick_params(axis='both', which='major', labelsize=24)
    ax.set_xlim(min_x_ang, max_x_ang)

    # Changing colorbar location setup
    divider = make_axes_locatable(ax)
    color_ax = divider.append_axes("bottom", size="5%", pad=1.5)
    color_bar = fig.colorbar(image, cax=color_ax, orientation='horizontal')
    color_bar.set_label('Counts', fontsize=30, labelpad=10)
    color_bar.ax.xaxis.set_label_position('bottom')
    color_bar.ax.xaxis.tick_top()
    color_bar.ax.tick_params(labelsize=24)  # Set the size of the colorbar tick labels

    plt.savefig(fig_path + f"2D_Spec_With_Noise_{target_name}_{num}.png", bbox_inches="tight")

    plt.show()


#########################################
# Plot reconstructed, low res SNR image #
#########################################
def plot_reconstructed_images(flipped_array, target_name):
    """
    Parameters
    ----------
    flipped_array: numpy array
        SNR brightness values in shape of image

    target_name: string
        specifies name of SN target

    Returns
    -------
    None: purpose is to plot reconstructed SNR image
    """

    #######################
    # Plot Original Image #
    #######################

    # Initialize figure
    fig, ax = plt.subplots(figsize=(7, 9))

    # Display the image using a magma colormap
    image = plt.imshow(flipped_array, cmap="viridis", aspect="auto")

    # Remove axis labels and ticks
    ax.set_title(f"{target_name} Reconstructed Image", fontsize=30)
    ax.set_xlabel("", fontsize=20)
    ax.set_ylabel("", fontsize=20)
    ax.tick_params(axis='both', which='both', bottom=False,
                   top=False, left=False, right=False, labelbottom=False, labelleft=False)

    # Changing colorbar location setup
    divider = make_axes_locatable(ax)
    color_ax = divider.append_axes("bottom", size="5%", pad=0.7)
    color_bar = fig.colorbar(image, cax=color_ax, orientation='horizontal')
    color_bar.set_label('Photon Counts', fontsize=30, labelpad=10)
    color_bar.ax.xaxis.set_label_position('bottom')
    color_bar.ax.xaxis.tick_top()
    color_bar.ax.tick_params(labelsize=24)  # Set the size of the colorbar tick labels

    plt.savefig(fig_path + f"{target_name}_reconstructed_at_1200_Ang.png", bbox_inches="tight")

    plt.show()
