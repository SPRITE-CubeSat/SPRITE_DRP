import matplotlib.pyplot as plt
import numpy as np
import os
import sys 
import argparse


def parse_arguments():
    """Handle argument input"""
    parser = argparse.ArgumentParser(description="Flux calibration using a standard star")

    # Required arguments
    parser.add_argument(
        "adcs_bin_directory",
        type=str,
        help="Directory containing the set of ADCS .bin output files to parse. Each bin file corresponds to an image\
        row. Make sure the directory only contains output for a single tracker image."
    )
    parser.add_argument(
        "--right_edge",
        action='store_true',
        help="Set this flag to manually denote that the right edge of the image is contained in the data. If so, the\
         last data-word on each line will only be parsed as containing two pixels. (This behaviour is automatic if the\
         line contains the full 427 bytes of a row.)",
    )
    return parser.parse_args()


def bytes_to_binary_string(byte_array, fmt='dw'):
    """Get more easily readable binary string for a byte array"""
    str_out = ''.join(format(byte, '08b')[::-1] for byte in byte_array[::-1])
    # Add spacing between pixels in dataword
    if fmt == 'dw':
        str_out = str_out[:10] + ' ' + str_out[10:20] + ' ' + str_out[20:30] + ' ' + str_out[30:]
    # Add spacing between each byte
    elif fmt == 'byte':
        str_out = str_out[:8] + ' ' + str_out[8:16] + ' ' + str_out[16:24] + ' ' + str_out[24:32]
    return str_out


def parse_adcs_dataword(dw_bytes):
    """Parse the bytes of a tracker image data word into pixel values"""
    # Get the binary string for the dataword
    bin_str = ''.join(format(byte, '08b') for byte in dw_bytes)

    # Third, extract the 0-9, 10-19, and 20-29 bits - but flip the order before casting to an int
    px1_b = bin_str[22:32]
    px2_b = bin_str[12:22]
    px3_b = bin_str[2:12]

    px1_int = int(px1_b, 2)
    px2_int = int(px2_b, 2)
    px3_int = int(px3_b, 2)

    return [px1_int, px2_int, px3_int]


def extract_tracker_image_from_bin_files(file_dir, right_edge=False):
    """Build up a tracker image from the .bin files saved by Hydra

    Args:
        file_dir (str): The absolute or relative path to the directory containing the .bin files saved by Hydra from
            an ADCS memory dump. Each .bin file contains an image row, which is 1708 bytes for a full row (1280px).
        right_edge (bool): Manual override to treat the right-most pixel as the right edge of the image. If so, the last
            data word in each line is parsed as only having two pixel values.

    Returns:
        numpy.ndarray: The tracker image
    """

    # Create a list structure to be built into a 2D tracker image by adding rows
    tracker_image = []

    # Iterate over the list of .bin files from Hydra, each corresponding to a singe 'readImageLine' command
    # These files are 2kB for a full image row
    for i, file_name in enumerate(sorted(os.listdir(file_dir))):
        with open(file_dir + '/' + file_name, 'rb') as f:

            # Get the bytes for this image row
            image_row_bytes = f.read()
            image_row_px = []

            # Make sure that all the image lines have the same number of bytes
            if i == 0:
                row_size_bytes = len(image_row_bytes)
            else:
                if len(image_row_bytes) != row_size_bytes:
                    raise RuntimeError("All .bin input files must contain the same number of bytes")

            # Iterate through the byte array 4 bytes at a time (32-bit "data word") and add the parsed pixel values
            for i in range(0, len(image_row_bytes), 4):
                image_row_px += parse_adcs_dataword(image_row_bytes[i:i + 4])

            # If we read a whole image line or the right-side of an image, we need to trim the last pixel
            if len(image_row_bytes) == 427 or right_edge:
                image_row_px = image_row_px[:-1]

            # Add the image row to the tracker image
            tracker_image.append(image_row_px)

    # Convert the tracker image from a 2D list to a 2D NumPy array
    tracker_image = np.array(tracker_image)

    return tracker_image


if __name__ == "__main__":

    args = parse_arguments()
    bin_files_dir = os.path.abspath(args.adcs_bin_directory)
    image_array = extract_tracker_image_from_bin_files(bin_files_dir)

    # Plot the image
    fig, axes = plt.subplots(1, 3, figsize=(16, 3), width_ratios=[5, 0.1, 1])
    im = axes[0].pcolor(image_array, cmap='gray', vmin=0, vmax=100)
    fig.suptitle(sys.argv[1])
    cb = fig.colorbar(im, cax=axes[1])
    axes[2].hist(image_array.flatten(), range=(0, np.max(image_array)), bins=int(np.max(image_array)))
    axes[2].set_yscale("log")
    cb.set_label("Counts")
    fig.tight_layout()
    fig.show()
    input("")