###############
# Directories #
###############

# path to files
file_path = "/Users/elenacarlson/PycharmProjects/Calibration_Stars/data_files"

# path to save figures
figure_path = "/Users/elenacarlson/PycharmProjects/Calibration_Stars/figures/"

###########################
# detector/slit constants #
###########################
arcsec_res = 100  # arcsec^2
exposure_time = 3000  # seconds

detector_size = 39   # mm
spectral_scale = 22  # Angstroms/mm

# slit location
slit_loc_ang = 1600  # Angstroms

# width of good data in fits file (in pixels)
col_pixels = 1380  # columns
row_pixels = 760  # rows

# mm per pixel for detector
mm_per_x_pixel = detector_size/col_pixels

# average aperture pixel location
ap_pixel_loc = 1062

# mm range before and after aperture
mm_after_ap = (col_pixels - ap_pixel_loc) * mm_per_x_pixel
mm_before_ap = ap_pixel_loc * mm_per_x_pixel

# convert to wavelength range before and after
ang_after_ap = mm_after_ap * spectral_scale
ang_before_ap = mm_before_ap * spectral_scale

# angstrom range
min_x_ang = slit_loc_ang - ang_before_ap
max_x_ang = slit_loc_ang + ang_after_ap
ang_range = max_x_ang - min_x_ang

# slit constants
x_pix = 40  # slit pixel width
y_pix = 760  # slit pixel length

##########################
# SNR and Star Constants #
##########################

# Point Source standard dev for gaussian
# 140 arcsecond/mm
arcsec_per_pix = mm_per_x_pixel * 140

# Invert for pixels per arcsecond
pix_per_arcsec = 1/arcsec_per_pix

# Multiply by angular size of 14 arcsecond
pix_std_dev = 14 * pix_per_arcsec

# Supernova Remnant size in arcsec
SN_width = 1.4 * 60
SN_height = 1.2 * 60

#########
# Noise #
#########

# Dark rate
'''
in a best case scenario:
dark_rate = 0.5 * 1e-4**2  #counts/micron^2/s
'''
dark_rate = 2 * 1e-4**2  # counts/micron^2/s

# noise per angstrom for a point source
point_source_pix_height = 150  # microns
angstrom_width = 45  # microns
noise = point_source_pix_height * angstrom_width * dark_rate  # counts/angstrom/second

###########################
# Dictionary of Constants #
###########################

# Constants used for convolution
convolve_constants = {max_x_ang: (slit_loc_ang + ang_after_ap), min_x_ang: (slit_loc_ang - ang_before_ap),
                      ang_range: (max_x_ang - min_x_ang), row_pixels: 760, col_pixels: 1380, SN_height: (1.2*60),
                      pix_std_dev: (14 * pix_per_arcsec)}
