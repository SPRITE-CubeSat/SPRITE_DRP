import argparse
import numpy as np 

from astropy.io import fits 
from matplotlib import pyplot as plt
from matplotlib import colors
from scipy.optimize import curve_fit, differential_evolution
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter
from scipy.signal import medfilt, find_peaks

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
        default=800,
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


def get_star_yrange(spec2d, sig_guess0=2.0, range_sig_width=3.0):
    """Get the ymin and ymax for the star's location along the y axis"""
    """"""
    yaxis = np.arange(spec2d.shape[0])
    yprof = np.sum(spec2d, axis=1)
    popt, pcov = curve_fit(
        gaussian, 
        yaxis,
        yprof, 
        p0 = [
            np.max(yprof),
            yaxis[np.nanargmax(yprof)],
            sig_guess0
        ]
    )
    ymin = int(round(popt[1] - 3 * popt[2]))
    ymax = int(round(popt[1] + 3 * popt[2]))
    return ymin, ymax

if __name__ == "__main__":
    
    # Handle user input
    args = parse_arguments()
    
    min_peak_width_px = 3 # pixels

    # Load the data
    obs_spec2d = fits.open(args.observed_2d_spectrum)[0].data 
    known_spec1d = fits.open(args.known_1d_spectrum)[1].data

    # Adjust the known spectrum to match SPRITE's resolution (assuming it's higher res.)
    sprite_resel_Angstrom = 1.3 
    sprite_resel_px = sprite_resel_Angstrom / (known_spec1d['WAVE'][1] - known_spec1d['WAVE'][0])
    
    known_spec1d_smooth = gaussian_filter(known_spec1d['FLUX'], 1.5 * sprite_resel_px / 2.355)
    known_spec1d_smooth -= medfilt(known_spec1d_smooth, 5001)
    known_spec1d_smooth /= np.max(known_spec1d_smooth)
    known_spec1d_interp_f = interp1d(known_spec1d['WAVE'], known_spec1d_smooth)

    # Get the location of the observed star and then collapse a 1D spectrum
    xaxis = np.arange(obs_spec2d.shape[1])
    ymin, ymax = get_star_yrange(obs_spec2d)
    obs_spec1d = np.sum(obs_spec2d[ymin:ymax, :], axis=0)
    obs_spec1d_backup = obs_spec1d.copy()
    obs_spec1d -= medfilt(obs_spec1d, kernel_size=101)
    obs_spec1d /= np.max(obs_spec1d)

    line_fit_thresh = 0.001

    # Example plot
    guess_wav = 925 + 0.65 * xaxis 
    known_spec1d_interp = known_spec1d_interp_f(guess_wav)
    fig, axes = plt.subplots(2, 1, figsize=(16, 8))
    axes[0].plot(guess_wav, obs_spec1d, 'k.-')
    axes[0].set_title("Observed spectrum (filtered)")
    axes[0].set_ylabel("Relative Counts")
    axes[0].fill_between(guess_wav, np.zeros_like(obs_spec1d), obs_spec1d > line_fit_thresh, alpha=0.25)
    axes[1].plot(guess_wav, known_spec1d_interp, 'k.-')
    axes[1].fill_between(guess_wav, np.zeros_like(known_spec1d_interp), known_spec1d_interp > line_fit_thresh, alpha=0.25)
    axes[1].set_title("Known spectrum (filtered / interpolated)")
    axes[1].set_ylabel("Relative Flux")
    axes[1].set_xlabel("Wavelength [Angstrom] (Initial guess solution)")
    for ax in axes:
        ax.set_yscale("log")
        ax.set_ylim([0.001, 1])
    fig.show()

    obs_spec1d_mask = obs_spec1d > line_fit_thresh
    res = 100
    wav0_min, wav0_max = 920, 960 
    dwav_min, dwav_max = 0.5, 0.7
    wav0_vals = np.linspace(wav0_min, wav0_max, res)
    dwav_vals = np.linspace(dwav_min, dwav_max, res)
    coeff_vals = np.zeros((res, res))
    for i, wav0 in enumerate(wav0_vals):
        for j, dwav in enumerate(dwav_vals):

            guess_wav = wav0 + xaxis * dwav
            try: known_spec1d_interp = interp1d(known_spec1d['WAVE'], known_spec1d_smooth)(guess_wav)
            except: continue
            
            corr_coeff = np.corrcoef(obs_spec1d_mask, known_spec1d_interp > line_fit_thresh)
            coeff_vals[i, j] = corr_coeff[0, 1]

    i_opt, j_opt = np.unravel_index(np.nanargmax(coeff_vals), coeff_vals.shape)
    wav0_opt = wav0_vals[i_opt]
    dwav_opt = dwav_vals[j_opt]
    wav_sol = wav0_opt + xaxis * dwav_opt
    known_spec1d_sol = interp1d(known_spec1d['WAVE'], known_spec1d['FLUX'])(wav_sol)
    print(i_opt, j_opt)
    print(wav0_opt, dwav_opt)


    fig, axes = plt.subplots(3, 1, figsize=(12, 12))
    axes[0].pcolor(obs_spec2d, vmax=10 * np.mean(obs_spec2d))
    axes[0].plot([0, obs_spec2d.shape[1]], [ymin, ymin],  'k--')
    axes[0].plot([0, obs_spec2d.shape[1]], [ymax, ymax],  'k--')
    axes[0].set_ylim([ymin - 100, ymax + 100])
    axes[1].pcolor(coeff_vals)
    xticks = [int(round(x)) for x in axes[1].get_xticks()[:-1]]
    yticks = [int(round(y)) for y in axes[1].get_yticks()[:-1]]
    axes[1].set_xticks(xticks, labels=["{0:.2f}".format(dwav_vals[x]) for x in xticks])
    axes[1].set_yticks(yticks, labels=["{0:.1f}".format(wav0_vals[y]) for y in yticks])
    axes[1].set_title("Cross-correlation Matrix")
    axes[1].set_xlabel("Angstroms per pixel")
    axes[1].set_ylabel("Pixel 0 Wavelength [Angstrom]")
    axes[2].plot(wav_sol, known_spec1d_sol * np.max(obs_spec1d_backup) / np.max(known_spec1d_sol), 'k.-', label="Known IUE/FUSE Spectrum")
    axes[2].set_yscale("log")
    axes[2].plot(wav_sol, obs_spec1d_backup, 'r--', label="Wavelength Solution")
    axes[2].legend()
    axes[2].set_ylim([1, np.max(obs_spec1d_backup)])
    fig.tight_layout()
    fig.show()
    input("")