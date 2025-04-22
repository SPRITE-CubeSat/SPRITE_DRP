###############
# Directories #
###############

# path to files
file_path = "/home/elca1097/LASP_SPRITE/data_files/"

# path to save figures
figure_path = "/home/elca1097/LASP_SPRITE/code/split_up_code/calibration/demo_figures/"

star_data_path = "/home/elca1097/LASP_SPRITE/standard_star_data/"

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
