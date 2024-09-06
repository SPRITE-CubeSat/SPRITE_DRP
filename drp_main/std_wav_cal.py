import argparse
import numpy as np 

from astropy.io import fits 
from matplotlib import pyplot as plt
from matplotlib import colors
from scipy.optimize import curve_fit, differential_evolution
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter
from scipy.signal import medfilt, find_peaks
from scipy.optimize import minimize

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
    parser.add_argument(
        "--x_center_min",
        type=float,
        default=None,
        help="Minimum value for dispersion center in x axis (default: midpoint - 20px)."
    )
    parser.add_argument(
        "--x_center_max",
        type=float,
        default=None,
        help="Maximum value for dispersion center in x axis (default: midpoint + 20px)."
    )
    parser.add_argument(
        "--w_center_min",
        type=float,
        default=1340.0,
        help="Minimum value for dispersion center in wavelength (default: 1340A)."
    )
    parser.add_argument(
        "--w_center_max",
        type=float,
        default=1390.0,
        help="Maximum value for dispersion center in wavelength (default:1390A)."
    )
    parser.add_argument(
        "--min_angstrom_per_px",
        type=float,
        default=0.5,
        help="Minimum value for linear spectral scale (angstrom per px, default 0.5)."
    )
    parser.add_argument(
        "--max_angstrom_per_px",
        type=float,
        default=0.7,
        help="Maximum value for linear spectral scale (angstrom per px, default 0.7)."
    )
    parser.add_argument(
        "--min_x2_coeff",
        type=float,
        default=0,
        help="Minimum coefficient of quadratic term of wavelength solution (Default: 0)"
    )
    parser.add_argument(
        "--max_x2_coeff",
        type=float,
        default=1e-3,
        help="Maximum coefficient of quadratic term of wavelength solution (Default: 1e-3)"
    )
    parser.add_argument(
        "--y_star_position",
        type=int,
        default=None,
        help="Initial guess for the location of the star's spectrum along the y axis (default: brightest pixel)"
    )
    parser.add_argument(
        "--y_star_sig_width",
        type=int,
        default=3.0,
        help="Width to sum stellar spectrum over in y-axis in terms of PSF sigma (Default: 3.0)"
    )
    parser.add_argument(
        "--med_filt_width",
        type=int,
        default=0.25,
        help="Width of median filter, as fraction of array size, for subtracting slowly varying continuum"
    )
    parser.add_argument(
        "--rel_peak_thresh",
        type=float,
        default=0.003,
        help="Threshold in terms of relative counts/flux to use a peak in cross-correlation (Default: 0.003)"
    )
    return parser.parse_args()


def gaussian(x, amplitude, mean, stddev):
    """Simple 1D gaussian"""
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))


def get_star_yrange(spec2d, y_guess=None, sig_guess=None, range_sig_width=3.0):
    """Get the ymin and ymax for the star's location along the y axis"""
    # Create y axis and 1D profile
    yaxis = np.arange(spec2d.shape[0])
    yprof = np.sum(spec2d, axis=1)

    # Get defaults
    if y_guess is None: y_guess = yaxis[np.nanargmax(yprof)]
    if sig_guess is None: sig_guess = 2.0 
    
    # Optimize values of gaussian fit
    popt, pcov = curve_fit(
        gaussian, 
        yaxis,
        yprof, 
        p0 = [
            np.max(yprof),
            y_guess,
            sig_guess
        ]
    )
    
    # Return the range in y values
    ymin = int(round(popt[1] - range_sig_width * popt[2]))
    ymax = int(round(popt[1] + range_sig_width * popt[2]))
    return ymin, ymax


def get_wavelength(x, params):
    """Quadratic solution"""
    x_rel = x - params[0]
    return params[1] * x_rel**2 + params[2] * x_rel + params[3]


def get_correlation_score(obs_spec_norm, known_spec_norm, rel_threshold, display=False):
    """Calculate correlation of two normalized spectra after thresholding"""
    obs_mask = obs_spec_norm >= rel_threshold
    known_mask = known_spec_norm >= rel_threshold
    corr_res = np.corrcoef(obs_mask, known_mask)[0, 1]

    
    if display:
        # Example plot
        fig, axes = plt.subplots(2, 1, figsize=(16, 8))
        fig.suptitle(corr_res)
        axes[0].plot(obs_spec_norm, 'k.-')
        axes[0].set_title("Observed spectrum (filtered)")
        axes[0].set_ylabel("Relative Counts")
        axes[0].fill_between(np.arange(len(obs_mask)), np.zeros_like(obs_mask), obs_mask, alpha=0.25)
        axes[1].plot(known_spec_norm, 'k.-')
        axes[1].fill_between(np.arange(len(known_mask)), np.zeros_like(known_mask), known_mask, alpha=0.25)
        axes[1].set_title("Known spectrum (filtered / interpolated)")
        axes[1].set_ylabel("Relative Flux")
        axes[1].set_xlabel("Wavelength [Angstrom] (Initial guess solution)")
        for ax in axes:
            ax.set_yscale("log")
            ax.set_ylim([0.001, 1])
        fig.show()
        input("")
        plt.close()

    return corr_res


def get_median_filter_width(data, fraction):
    """Get an integer, odd median filter width"""
    width = int(len(data * fraction) // 2.0)
    if width % 2 == 0: return width + 1
    return width

# Example objective function
def objective_function(params, x_observed, counts_observed, wav_known, flux_known, rel_thresh):
    """Objective function for wavelenght solution fitting: negative correlation with known spec"""
    # 1. Get the wavelength solution using the current parameters
    wav_sol = get_wavelength(x_observed, params)
    
    # 2. Interpolate the known spectrum onto the new wavelength grid
    interp_flux = np.interp(wav_sol, wav_known, flux_known)  # Interpolation
    
    # 3. Compute the correlation score
    xcor = get_correlation_score(counts_observed, interp_flux, rel_threshold=rel_thresh)

    # Since most optimizers minimize functions, we return the negative of the correlation score
    return -xcor

# Optimize the parameters using differential_evolution
def optimize_wavelength_solution_de(x_observed, counts_observed, wav_known, flux_known, param_bounds, rel_thresh=0.003):
    """
    Optimizes the wavelength solution parameters using differential evolution.
    
    Parameters:
    x_observed (array-like): The x-axis of the observed spectrum.
    counts_observed (array-like): The observed spectrum counts.
    wav_known (array-like): The known wavelength grid from the reference spectrum.
    flux_known (array-like): The known flux values from the reference spectrum.
    param_bounds (list of tuple): Bounds for each parameter in the form [(min1, max1), (min2, max2), ...].
    
    Returns:
    result (OptimizeResult): The result of the optimization process.
    """
    # Use scipy.optimize.differential_evolution to optimize the params
    result = differential_evolution(objective_function, bounds=param_bounds, 
                                    args=(x_observed, counts_observed, wav_known, flux_known, rel_thresh))
    
    return result

if __name__ == "__main__":
    
    # Handle user input
    args = parse_arguments()

    # Load the data
    obs_spec2d = fits.open(args.observed_2d_spectrum)[0].data 
    known_spec1d = fits.open(args.known_1d_spectrum)[1].data

    # Smooth the known spectrum to match SPRITE's resolution (assuming it's higher res.)
    print("Preparing known spectrum... (TODO: Speed this up by pre-processing the spectrum)")
    sprite_resel_Angstrom = 1.3 
    known_wav = known_spec1d['WAVE']
    sprite_resel_px = sprite_resel_Angstrom / (known_wav[1] - known_wav[0])
    known_spec1d_smooth = gaussian_filter(known_spec1d['FLUX'], 5 * sprite_resel_px / 2.355)
    
    # Subtract background (median fitler) and then normalize the spectrum
    known_spec1d_bg =  medfilt(known_spec1d_smooth, 7001)# get_median_filter_width(known_spec1d_smooth, args.med_filt_width))
    known_spec1d_bgsub = known_spec1d_smooth - known_spec1d_bg
    known_spec1d_norm = known_spec1d_bgsub / np.max(known_spec1d_bgsub)

    # Get the location of the observed star and then collapse a 1D spectrum
    print("Extracting observed spectrum...")
    xaxis = np.arange(obs_spec2d.shape[1])
    ymin, ymax = get_star_yrange(obs_spec2d)
    obs_spec1d = np.sum(obs_spec2d[ymin:ymax, :], axis=0)

    # Filter the spectrum
    obs_spec1d_bg = medfilt(obs_spec1d, 401)# get_median_filter_width(obs_spec1d, args.med_filt_width))
    obs_spec1d_bgsub = obs_spec1d - obs_spec1d_bg
    obs_spec1d_norm = obs_spec1d_bgsub / np.max(obs_spec1d_bgsub)

    # Create bounds on each parameter for differential evolution 
    param_bounds = [
        # Reference x-pixel at which center of dispersion lies
        (
            len(xaxis)/2.0 - 30 if args.x_center_min is None else args.x_center_min,
            len(xaxis)/2.0 + 30 if args.x_center_min is None else args.x_center_max
        ),
        # Coefficient of x2 term in the wavelength solution
        (
            args.min_x2_coeff, 
            args.max_x2_coeff
        ),
        # Linear angstrom per pixel term
        (
            args.min_angstrom_per_px,
            args.max_angstrom_per_px
        ),
        # Refernece wavelength at center of dispersion
        (
            args.w_center_min,
            args.w_center_max
        )
    ]
  
    # Run optimization
    result = optimize_wavelength_solution_de(xaxis, obs_spec1d_norm, known_wav, known_spec1d_norm, param_bounds, rel_thresh=args.rel_peak_thresh)

    # Optimized parameters
    optimized_params = result.x
    optimized_wav = get_wavelength(xaxis, optimized_params)
    optimized_spec = np.interp(optimized_wav, known_wav, known_spec1d_norm)
    print("Optimized parameters:", optimized_params)

    # Print out table of solved values
    print("{0:5}\t{1:10}".format("x", "wav_Angstrom"))
    for i, wav in enumerate(optimized_wav):
        print(f"{i:5}\t{wav:10.2f}")

    # Print out solution equation
    xref, c2, c1, c0 = optimized_params
    print("wav(x) = {0:.2f} + {1:.2E} * (x - {2:.2f})^2 + {3:.2f} * (x - {2:.2f})".format(c0, c2, xref, c1))


    # Plot the result
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes[0, 0].pcolor(obs_spec2d, vmax=10 * np.mean(obs_spec2d))
    axes[0, 0].plot([0, obs_spec2d.shape[1]], [ymin, ymin],  'k--')
    axes[0, 0].plot([0, obs_spec2d.shape[1]], [ymax, ymax],  'k--')
    axes[0, 0].set_ylim([ymin - 100, ymax + 100])
    axes[1, 0].plot(known_wav, known_spec1d_smooth * np.max(obs_spec1d) / np.max(known_spec1d_smooth), 'k.-', label="Known IUE/FUSE Spectrum")
    axes[1, 0].set_yscale("log")
    axes[1, 0].plot(optimized_wav, obs_spec1d, 'r-', label="Wavelength Solution")
    axes[1, 0].legend()
    axes[1, 0].set_ylim([1, np.max(obs_spec1d)])
    axes[1, 0].plot(optimized_wav, obs_spec1d_bg, 'r:')
    axes[1, 0].plot(known_wav, known_spec1d_bg * np.max(obs_spec1d) / np.max(known_spec1d_smooth), 'k:')
    axes[0, 1].plot(known_wav, known_spec1d_norm, 'k.-', label="Normalized/filtered Known Spec")
    axes[0, 1].plot([known_wav[0], known_wav[-1]], [args.rel_peak_thresh]*2, 'r--')
    axes[0, 1].fill_between(known_wav, np.zeros_like(known_wav), known_spec1d_norm > args.rel_peak_thresh, alpha=0.25)
    axes[1, 1].plot(optimized_wav, obs_spec1d_norm, 'k.-', label="Normalized/filtered Observed Spec")
    axes[1, 1].plot([optimized_wav[0], optimized_wav[-1]], [args.rel_peak_thresh]*2, 'r--')
    axes[1, 1].fill_between(optimized_wav, np.zeros_like(optimized_wav), obs_spec1d_norm > args.rel_peak_thresh, alpha=0.25)
    axes[1, 1].sharex(axes[0, 1])
    for ax in axes[:, 1]:
        #ax.set_xlim([900, 2000])
        ax.legend()
        ax.set_yscale("log")
        ax.set_ylim([args.rel_peak_thresh / 2.0, 1])
    fig.tight_layout()
    fig.show()
    input("")