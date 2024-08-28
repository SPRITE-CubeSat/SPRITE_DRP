import argparse
import numpy as np 

from astropy.io import fits 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

def parse_arguments():
    parser = argparse.ArgumentParser(description="Wavelength calibration using a standard star")

    # Required arguments
    parser.add_argument(
        "observed_2d_spectrum",
        type=str,
        help="Path to the FITS file containing the observed 2D spectrum of the standard star"
    )
    parser.add_argument(
        "known_1d_spectrum",
        type=str,
        help="Path to the FITS file containing the known 1D spectrum of the standard star"
    )

    # Optional arguments with defaults
    parser.add_argument(
        "--min_wavelength",
        type=float,
        default=950,
        help="Minimum wavelength in Angstroms to use for calibration (default: 900)"
    )
    parser.add_argument(
        "--angstrom_per_pixel",
        type=float,
        default=0.6,
        help="Initial guess for the angstrom per pixel scale of the image (default: 0.2)"
    )
    parser.add_argument(
        "--y_position",
        type=int,
        default=1024,
        help="Initial guess for the location of the star's spectrum along the y axis (default: 1024)"
    )

    return parser.parse_args()


def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))


if __name__ == "__main__":
    args = parse_arguments()
    
    obs_spec2d = fits.open(args.observed_2d_spectrum)[0].data 
    known_spec1d = fits.open(args.known_1d_spectrum)[1].data

    yaxis = np.arange(obs_spec2d.shape[0])
    yprof = np.sum(obs_spec2d, axis=1)

    popt, pcov = curve_fit(
        gaussian, 
        yaxis,
        yprof, 
        p0 = [
            np.max(yprof),
            yaxis[np.nanargmax(yprof)],
            2.0
        ]
    )
    ymin = int(round(popt[1] - 3 * popt[2]))
    ymax = int(round(popt[1] + 3 * popt[2]))
    obs_spec1d = np.sum(obs_spec2d[ymin:ymax, :], axis=0)

    xaxis = np.arange(obs_spec2d.shape[1])
    guess_wav = xaxis * args.angstrom_per_pixel + args.min_wavelength 

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes[0, 0].pcolor(obs_spec2d)
    axes[1, 0].plot(known_spec1d["WAVE"], known_spec1d["FLUX"], 'k.-', label="Known IUE/FUSE Spectrum")
    axes[1, 0].set_yscale("log")
    axes_1_0_twinx = axes[1, 0].twinx()
    axes_1_0_twinx.plot(guess_wav, obs_spec1d, 'r--', label="Initial Guess Solution")
    axes_1_0_twinx.legend()
    axes_1_0_twinx.set_yscale("log")
    axes[1, 0].legend(loc=5)
    axes[0, 1].plot(yprof, yaxis, 'k.-', label="Y Profile")
    axes[0, 1].plot(gaussian(yaxis, *popt), yaxis, 'r-', label="Gaussian Fit")
    axes[0, 1].legend()

    fig.show()
    input("")