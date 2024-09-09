from astropy.io import fits
from scipy.interpolate import interp1d 

import matplotlib.pyplot as plt 
import numpy as np 

import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Flux calibration using a standard star")

    # Required arguments
    parser.add_argument(
        "fuse_spectrum",
        type=str,
        help="Path to the FITS file containing the FUSE spectrum for the standard star"
    )
    parser.add_argument(
        "iue_spectrum",
        type=str,
        help="Path to the FITS file containing the IUE spectrum from the standard star"
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="What to save the output file as"
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    # Load and extract data
    fuse_data = fits.open(args.fuse_spectrum)[1].data
    iue_data = fits.open(args.iue_spectrum)[1].data
    fuse_wav, fuse_flux = fuse_data["WAVE"][0], fuse_data["FLUX"][0]
    iue_wav, iue_flux = iue_data["WAVE"][0], iue_data["FLUX"][0]

    # Get linear interpolation of each
    iue_interp_f = interp1d(iue_wav, iue_flux, kind='linear')
    fuse_interp_f = interp1d(fuse_wav, fuse_flux, kind='linear')

    # Create full wav axis with 
    merge_dlam = fuse_wav[1] - fuse_wav[0]
    merge_wav = np.arange(fuse_wav[0], iue_wav[-1], step=0.5)
    iue_indices = merge_wav >= iue_wav[0]
    fuse_indices = merge_wav <= fuse_wav[-1]
    overlap_indices = (iue_indices & fuse_indices)

    # Create interpolated spectra
    fuse_interp = np.zeros_like(merge_wav)
    fuse_interp[fuse_indices] = fuse_interp_f(merge_wav[fuse_indices])
    iue_interp = np.zeros_like(merge_wav)
    iue_interp[iue_indices] = iue_interp_f(merge_wav[iue_indices])

    # Create the merged spectrum
    merge_spec = np.zeros_like(merge_wav)
    merge_spec[fuse_indices] = fuse_interp[fuse_indices]
    merge_spec[iue_indices] = iue_interp[iue_indices]
    merge_spec[overlap_indices] = (fuse_interp[overlap_indices] + iue_interp[overlap_indices]) / 2.0

    # Save to FITS
    spec1d_hdu = fits.BinTableHDU.from_columns(
        fits.ColDefs([
            fits.Column(
                name="WAVE",
                format='D', #16-bit Integer = "I",
                unit="Angstrom",
                array = merge_wav
            ),
            fits.Column(
                name="FLUX",
                format='D', #16-bit Integer = "I",
                unit="erg/cm^2/s/Angstrom",
                array = merge_spec
            )
        ])
    )

    spec1d_fits = fits.HDUList([fits.PrimaryHDU(), spec1d_hdu])
    spec1d_fits.writeto(args.output_file, overwrite=True, checksum=True)

    fig, ax = plt.subplots(1, 1)
    ax.plot(merge_wav, iue_interp, label="IUE")
    ax.plot(merge_wav, fuse_interp, label="FUSE")
    ax.plot(merge_wav, merge_spec, label="MERGED")
    ax.set_xlabel("Wavelength [Angstrom]")
    ax.set_ylabel("Flux [erg/s/cm2/Angstrom]")
    ax.legend()
    fig.show()
    input("")
