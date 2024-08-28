from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

import numpy as np
import time
import meta_data
import pathlib

# Always save files to sample_data/generated/ regardless of where the script is run
generated_dir = str(pathlib.Path(__file__).parent.resolve()) + "/generated/"

#################
# PHOTLIST.FITS #
#################
photlist_filename = meta_data.get_mast_filename(setting='ttag', filetype='photlist.fits')

# Params for the dummy data
table_size = 10  # Number of photon events

# Create columns
axis1_col = fits.Column(
    name="X",
    format='I',  # 16-bit Integer = "I",
    unit="pix",
    array=np.random.randint(low=1, high=2048, size=table_size, dtype=int)
)
axis2_col = fits.Column(
    name="Y",
    format='I',  # 16-bit Integer = "I",
    unit="pix",
    array=np.random.randint(low=1, high=2048, size=table_size, dtype=int)
)
ph_col = fits.Column(
    name="PULSE_HEIGHT",
    format='B',  # 8-bit Integer = "B",
    unit="adu",
    array=np.random.randint(low=1, high=256, size=table_size, dtype=int)
)
time_col = fits.Column(
    name="TIME",
    format='E',  # Single-precision (32-bit) float = "E"
    unit='s',
    array=np.random.rand(table_size) * meta_data.T_EXP
)
cols = fits.ColDefs([axis1_col, axis2_col, ph_col, time_col])

# Create primary HDU and table HDU
photlist_hdu = fits.BinTableHDU.from_columns(cols)
primary_hdu = fits.PrimaryHDU()

# Add keywords from certain required collections
phdu_dicts_to_add = [
    meta_data.common_mast_header_dict,
    meta_data.spec2d_mast_header_dict,
    meta_data.sprite_housekeeping_dict
]
for dict_ in phdu_dicts_to_add:
    for key, val in dict_.items():
        if key not in primary_hdu.header.keys():
            primary_hdu.header[key] = val

# Update obs mode to TTAG
primary_hdu.header["OBS_MODE"] = "TTAG"
primary_hdu.header["DAT_TYPE"] = "PHOTLIST"

# Add helpful comments for each TTYPE field
for i, (key, val) in enumerate(photlist_hdu.header.items()):
    if key[:5] != "TTYPE": continue
    if val == "X":
        photlist_hdu.header[key] = (val, "Horizontal (spectral) axis (AXIS1)")
    elif val == "Y":
        photlist_hdu.header[key] = (val, "Vertical (spatial) axis (AXIS2)")
    elif val == "PULSE_HEIGHT":
        photlist_hdu.header[key] = (val, "Pulse height of photon event on MCP")
    elif val == "TIME_DELTA":
        photlist_hdu.header[key] = (val, "Time (s) since MJD-OBS in table Header")

# Add MJDREF column to Table header for time column
photlist_hdu.header["JDREF"] = (
    primary_hdu.header["MJD-BEG"], "Reference MJD for TIME column")

# Create and export FITS
photlist_fits = fits.HDUList([primary_hdu, photlist_hdu])
photlist_fits.writeto(generated_dir + photlist_filename,
                      overwrite=True, checksum=True)
print(f"Saved SPRITE_DRP/sample_data/generated/{photlist_filename}")

##############
# INT2D.FITS #
##############
int2d_filename = meta_data.get_mast_filename(setting='accum', filetype='int2d.fits')

# Create the dummy 2D spec image data
int2d_data = np.zeros((meta_data.X_SIZE, meta_data.Y_SIZE)).astype(np.int32)

# Create the HDU
int2d_hdu = fits.PrimaryHDU(int2d_data)

# Add keywords from certain required collections
dicts_to_include = [
    meta_data.common_mast_header_dict,
    meta_data.spec2d_mast_header_dict,
    meta_data.sprite_housekeeping_dict
]
for dict_ in dicts_to_include:
    for key, val in dict_.items():
        if key not in int2d_hdu.header.keys():
            int2d_hdu.header[key] = val

# Update obs mode to TTAG
int2d_hdu.header["OBS_MODE"] = "ACCUM"

# Create and write the FITS file
int2d_fits = fits.HDUList([int2d_hdu])
int2d_fits.writeto(generated_dir + int2d_filename, overwrite=True, checksum=True)
print(f"Saved SPRITE_DRP/sample_data/generated/{int2d_filename}")

###############
# PHDIST.FITS #
###############
# @title _phdist.fits {vertical-output:true}
phdist_filename = meta_data.get_mast_filename(setting='accum', filetype='phdist.fits')

# Simulate PH data
ph_vals = np.arange(0, 256)
ph_counts = (1000 * np.exp(
    -(ph_vals - 150) ** 2 / (2 * 5 ** 2)
) + np.random.poisson(1, 256)
             ).astype(int)

# Create columns
axis1_col = fits.Column(
    name="PH",
    format='B',  # 8-bit Integer = "B",
    unit="adu",
    array=ph_vals
)
axis2_col = fits.Column(
    name="Count",
    format='I',  # 16-bit Integer = "I",
    unit="photon",
    array=ph_counts
)
cols = fits.ColDefs([axis1_col, axis2_col])

# Create primary HDU and table HDU
phdist_hdu = fits.BinTableHDU.from_columns(cols)
primary_hdu = fits.PrimaryHDU()

# Add keywords from certain required collections
phdu_dicts_to_add = [
    meta_data.common_mast_header_dict,
    meta_data.sprite_housekeeping_dict
]
for dict_ in phdu_dicts_to_add:
    for key, val in dict_.items():
        if key not in primary_hdu.header.keys():
            primary_hdu.header[key] = val

# Create and export FITS
phdist_fits = fits.HDUList([primary_hdu, phdist_hdu])
phdist_fits.writeto(generated_dir + phdist_filename, overwrite=True, checksum=True)
print(f"Saved SPRITE_DRP/sample_data/generated/{phdist_filename}")

###############
# CALIMG.FITS #
###############
calimg_filename = meta_data.get_mast_filename(setting='accum', filetype='calimg.fits')
calimg_data = np.zeros((meta_data.CALIMG_X_SIZE, meta_data.CALIMG_Y_SIZE))
calimg_hdu = fits.PrimaryHDU(calimg_data)

# Add keywords from certain required collections
dicts_to_include = [
    meta_data.common_mast_header_dict,
    meta_data.image_mast_header_dict,
    meta_data.sprite_housekeeping_dict
]
for dict_ in dicts_to_include:
    for key, val in dict_.items():
        if key not in calimg_hdu.header.keys():
            calimg_hdu.header[key] = val

# Update obs mode to ACCUM/CAL
calimg_hdu.header["OBS_MODE"] = "ACCUM"
calimg_hdu.header["DAT_TYPE"] = "CALIMG"
# Create and write the FITS file
calimg_fits = fits.HDUList([calimg_hdu])
calimg_fits.writeto(generated_dir + calimg_filename, overwrite=True, checksum=True)
print(f"Saved SPRITE_DRP/sample_data/generated/{calimg_filename}")

###############
# SPEC1D.FITS #
###############
spec1d_filename = meta_data.get_mast_filename(setting='accum', filetype='spec1d.fits')

# Create columns
axis1_col = fits.Column(
    name="WAVE",
    format='I',  # 8-bit Integer = "B",
    unit="Angstrom",
    array=np.linspace(900, 1900, 2048)
)
axis2_col = fits.Column(
    name="FLUX",
    format='I',  # 16-bit Integer = "I",
    unit="erg/cm^2/s/Angstrom",
    array=np.zeros(2048)
)
axis3_col = fits.Column(
    name="FLUX_ERR",
    format='I',  # 16-bit Integer = "I",
    unit="erg/cm^2/s/Angstrom",
    array=np.zeros(2048)
)
cols = fits.ColDefs([axis1_col, axis2_col, axis3_col])

# Create primary HDU and table HDU
spec1d_hdu = fits.BinTableHDU.from_columns(cols)
primary_hdu = fits.PrimaryHDU()

# Add keywords from certain required collections
dicts_to_include = [
    meta_data.common_mast_header_dict,
    meta_data.sprite_housekeeping_dict
]
for dict_ in dicts_to_include:
    for key, val in dict_.items():
        if key not in primary_hdu.header.keys():
            primary_hdu.header[key] = val

# Create and write the FITS file
spec1d_fits = fits.HDUList([primary_hdu, spec1d_hdu])
spec1d_fits.writeto(generated_dir + spec1d_filename, overwrite=True, checksum=True)
print(f"Saved SPRITE_DRP/sample_data/generated/{spec1d_filename}")

###############
# SPEC2D.FITS #
###############
spec2d_filename = meta_data.get_mast_filename(setting='accum', filetype='spec2d.fits')

# Create the dummy 2D spec image data
spec2d_data = np.zeros((meta_data.X_SIZE, meta_data.Y_SIZE)).astype(np.int32)

# Create the HDU
primary_hdu = fits.PrimaryHDU()
spec2d_hdu = fits.ImageHDU(spec2d_data)
spec2d_err_hdu = fits.ImageHDU(spec2d_data)

# Add keywords from certain required collections
dicts_to_include = [
    meta_data.common_mast_header_dict,
    meta_data.spec2d_mast_header_dict,
    meta_data.sprite_housekeeping_dict
]
for dict_ in dicts_to_include:
    for key, val in dict_.items():
        if key not in spec2d_hdu.header.keys():
            primary_hdu.header[key] = val

# Name the extensions
primary_hdu.header["EXTNAME"] = "Primary"
spec2d_hdu.header["EXTNAME"] = "SCI"
spec2d_err_hdu.header["EXTNAME"] = "ERR"

# Add brightness units
spec2d_hdu.header["OBS_MODE"] = "ACCUM"
spec2d_hdu.header["BUNIT"] = "erg/s/cm^2/arcsec^2/Angstrom"
spec2d_err_hdu.header["OBS_MODE"] = "ACCUM"
spec2d_err_hdu.header["BUNIT"] = "erg/s/cm^2/arcsec^2/Angstrom"

# Create and write the FITS file
spec2d_fits = fits.HDUList([primary_hdu, spec2d_hdu, spec2d_err_hdu])
spec2d_fits.writeto(generated_dir + spec2d_filename, overwrite=True, checksum=True)
print(f"Saved SPRITE_DRP/sample_data/generated/{spec2d_filename}")