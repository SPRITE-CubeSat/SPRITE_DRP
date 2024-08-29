###############
# Directories #
###############

# path to files
file_path = "/Users/elenacarlson/PycharmProjects/SNR_files/files/"

# path to save figures
fig_path = "/Users/elenacarlson/Desktop/"

###########################
# detector/slit constants #
###########################

detector_size = 39    # mm
spectral_scale = 22   # Angstroms/mm

# Slit location
slit_loc_ang = 1600   # Angstroms

# Width of good data in fits file (in pixels)
col_pixels = 1380   # columns
row_pixels = 760    # rows

# mm per pixel for detector
mm_per_x_pixel = detector_size/col_pixels

# Average aperture pixel location
# Determined prior using for loop
ap_pixel_loc = 1062

# mm range before and after aperture
mm_after_ap = (col_pixels - ap_pixel_loc) * mm_per_x_pixel
mm_before_ap = ap_pixel_loc * mm_per_x_pixel

# Convert to wavelength range before and after
ang_after_ap = mm_after_ap * spectral_scale
ang_before_ap = mm_before_ap * spectral_scale

# Angstrom range
min_x_ang = slit_loc_ang - ang_before_ap
max_x_ang = slit_loc_ang + ang_after_ap
ang_range = max_x_ang - min_x_ang

# Slit constants
x_pix = 40         # Slit pixel width
y_pix = 760        # Slit pixel length

# Quantum Efficiency
quantum_efficiency = 0.5

###################
# Other Constants #
###################

# Exposure Time in seconds
exposure_time = 3000

# Point Source standard dev for gaussian
# 140 arc-seconds / mm
arcsec_per_pix = mm_per_x_pixel * 140

# Invert for pixels per arc-seconds
pix_per_arcsec = 1/arcsec_per_pix

# Multiply by angular size of 14 arc-seconds
pix_std_dev = 14 * pix_per_arcsec

# Dark rate
dark_rate = 3        # Photon counts/sec

# Supernova Remnant size in arc-seconds
supernova_width = 1.4 * 60
supernova_height = 1.2 * 60

###########################
# Dictionary of Constants #
###########################

# Constants used to convolve data
convolve_constants = {max_x_ang: (slit_loc_ang + ang_after_ap), min_x_ang: (slit_loc_ang - ang_before_ap),
                      ang_range: (max_x_ang - min_x_ang), row_pixels: 760, col_pixels: 1380, supernova_height: (1.2*60),
                      pix_std_dev: (14 * pix_per_arcsec)}
