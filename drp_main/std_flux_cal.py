import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Flux calibration using a standard star")

    # Required arguments
    parser.add_argument(
        "observed_spectrum",
        type=str,
        help="Path to the FITS file containing the observed 2D spectrum of the standard star"
    )
    parser.add_argument(
        "known_spectrum",
        type=str,
        help="Path to the FITS file containing the known 1D spectrum of the standard star"
    )

    # Optional arguments with defaults
    parser.add_argument(
        "--min_wavelength",
        type=float,
        default=900,
        help="Minimum wavelength in Angstroms to use for calibration (default: 900)"
    )
    parser.add_argument(
        "--max_wavelength",
        type=float,
        default=1900,
        help="Maximum wavelength in Angstroms to use for calibration (default: 1900)"
    )
    parser.add_argument(
        "--angstrom_per_pixel",
        type=float,
        default=0.2,
        help="Initial guess for the angstrom per pixel scale of the image (default: 0.2)"
    )
    parser.add_argument(
        "--y_position",
        type=int,
        default=1024,
        help="Initial guess for the location of the star's spectrum along the y axis (default: 1024)"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    print(args)