from datetime import datetime, timedelta
from astropy.time import Time


def get_mast_filename(mission='sprite', setting='accum', target_name='m81',
                      obs_sequence=1, obs_number=1, wave='fuv', filetype='spec2d.fits',
                      version='1'):
    """Return a MAST compliant filename based on various parameters"""
    target_string = f"{target_name}-s{obs_sequence:02d}-n{obs_number:05d}"
    return f"mccm_{mission}_{setting}_{target_string}_{wave}_{version}_{filetype}"


# Parameters for the mock data
RA = 149.688833333
DEC = 47.0569166667
T_EXP = 100  # seconds
MIN_WAV = 900  # Angstrom
MAX_WAV = 1900  # Angstrom
X_SIZE = 2048
Y_SIZE = 2048
CALIMG_X_SIZE = 300
CALIMG_Y_SIZE = 300
#
#  Create datetime objects for the main header fields
#
datetime_begin = datetime.now()  # Time("2023-12-12 00:00:00", scale="utc")
datetime_mid = datetime_begin + timedelta(seconds=T_EXP / 2.0)
datetime_end = datetime_begin + timedelta(seconds=T_EXP)

# ISO Format String Representations of the dates
date_beg_iso = datetime_begin.isoformat(timespec='seconds')
date_mid_iso = datetime_end.isoformat(timespec='seconds')
date_end_iso = datetime_end.isoformat(timespec='seconds')

# MJD Representations of the dates
date_beg_mjd = Time(date_beg_iso, format="isot", scale="utc").mjd
date_mid_mjd = Time(date_mid_iso, format="isot", scale="utc").mjd
date_end_mjd = Time(date_end_iso, format="isot", scale="utc").mjd

# MAST Common Fields
common_mast_header_dict = {
    "DATE": (date_beg_iso, "File creation date"),
    "DATE-BEG": (date_beg_iso, "Observation start datetime"),
    "DATE-END": (date_end_iso, "Observation end datetime"),
    "TIMESYS": ("UTC", "Time system for ISO Dates"),
    "DOI": ("doi:10.17909/gv67-bx83", "SPRITE's MAST Homepage"),
    "TARG": ("M81", "Target name"),
    "VER": ("v1.0", "File version"),
    "INSTRUME": ("SPRITE", "Instrument"),
    "LICENSE": ("CC BY 4.0", ""),
    "LICENURL": ("https://creativecommons.org/licenses/by/4.0/", ""),
    "MJD-BEG": (date_beg_mjd, "Observation start as MJD"),
    "MJD-END": (date_end_mjd, "Observation end as MJD"),
    "MJD-MID": (date_mid_mjd, "Observation midpoint as MJD"),
    "OBSERVAT": ("SPRITE", "Observatory - always 'SPRITE' for SPRITE"),
    "TELAPSE": (T_EXP, "Elapsed time between obervation start and end"),
    "TELESCOP": ("SPRITE", "Telescope"),
    "INTENT": ("SCIENCE", "SCIENCE or CALIBRATION"),
    "XPOSURE": (T_EXP, "Exposure time minus dead time and lost time"),
    "RA_TARG": (RA, "Target (pointing) coordinates, ICRS deg"),
    "DEC_TARG": (DEC, "Target (pointing) coordinates, ICRS deg"),
    "WAVEBAND": ("UV", "Waveband in MAST DB"),
    "WAVE_MIN": (MIN_WAV, "Minimum wavelength in Angstrom"),
    "WAVE_MAX": (MAX_WAV, "Maximum wavelength in Angstrom"),
}

# MAST 2D Spec Fields
spec2d_mast_header_dict = {
    "BUNIT": ("photon", "Brightness unit"),
    "WCSAXES": (3, "Number of WCS Axes"),
    "RADESYS": ("ICRS", "RA-DEC System"),
    "FILTER": ("NONE", "Required MAST keyword but no filters in SPRITE"),
    # WAVE is WCS axis 1
    "CTYPE1": ("WAVE", "Wavelength"),
    "CUNIT1": ("Angstrom", "Wavelength units"),
    "CRVAL1": (1400.0, "WAVE Reference value (angstrom)"),
    # RA is WCS axis 2
    "CTYPE2": ("RA--TAN", "Right-ascension"),
    "CUNIT2": ("deg", "Right-ascension units"),
    "CRVAL2": (RA, "RA reference value (degrees)"),
    # DEC is WCS axis 3
    "CTYPE3": ("DEC--TAN", "Declination"),
    "CUNIT3": ("deg", "Declination units"),
    "CRVAL3": (DEC, "DEC reference value (degrees)"),
    # Reference pixels
    "CRPIX1": (1, "Reference X pixel"),
    "CRPIX2": (1, "Reference Y pixel"),
    "CRPIX3": (1, "Degenerate j=3 WCS pixel axis"),
    # Only WAVE changes along X axis
    "CD1_1": (0.2, "WAVE Angstrom per X pixel"),
    "CD2_1": (0, "RA degrees per X pixel"),
    "CD3_1": (0, "DEC degrees per X pixel"),
    # Only RA/DEC change along Y axis
    "CD1_2": (0, "WAVE Angstrom per Y pixel"),
    "CD2_2": (0.1 / 3600, "RA degrees per Y pixel"),
    "CD3_2": (-0.1 / 3600, "DEC degrees per Y pixel"),
    # Inert values required by FITS Standard
    "CD2_3": (1, "Inert value required by FITS Standard"),
    "CD3_3": (1, "Inert value required by FITS Standard")
}

# MAST 2D Image Fields
image_mast_header_dict = {
    "BUNIT": ("photon", "Brightness unit"),
    "WCSAXES": (2, "Number of WCS Axes"),
    "RADESYS": ("ICRS", "RA-DEC System"),
    "CTYPE1": ("RA---TAN", "Right-Ascension"),
    "CUNIT1": ("deg", "Degrees"),
    "CTYPE2": ("DEC--TAN", "Declination"),
    "CUNIT2": ("deg", "Degrees"),
    "CRVAL1": (RA, "RA reference value (degrees)"),
    "CRVAL2": (DEC, "DEC reference value (degrees)"),
    "CRPIX1": (1, "Reference X pixel"),
    "CRPIX2": (1, "Reference Y pixel"),
    "CD1_1": (8.243530717e-12, "RA degrees per X pixel"),
    "CD2_1": (0.00018859, "DEC degrees per X pixel"),
    "CD1_2": (8.096e-05, "RA degrees per Y pixel"),
    "CD2_2": (-3.53887399569e-12, "DEC degrees per Y pixel"),
    "FILTER": ("NONE", "Required MAST keyword but no filters in SPRITE")
}

# SPRITE House-keeping Fields
sprite_housekeeping_dict = {
    "SPRITEID": (1, "Unique SPRITE observation ID"),
    "STD_ID": (0, "ObsID of Standard Star to use for this file"),
    "ERR_CODE": ("", "SPRITE Flight-side Error Codes"),
    "DAT_TYPE": ("SPEC2D", "SPRITE Data Type"),
    "OBS_MODE": ("", "SPRITE Observing Mode"),
    "HV_T_AVG": (0.0, "High-voltage avg temp (degrees celsius)"),
    "HV_T_VAR": (0.0, "High-voltafe temp var (degrees celsius)"),
    "SM_T_AVG": (0.0, "Secondary avg temp (degrees celsius)"),
    "SM_T_VAR": (0.0, "Secondary temp var (degrees celsius)"),
    "DT_T_AVG": (0.0, "Detector avg temp (degrees celsius)"),
    "DT_T_VAR": (0.0, "Detector temp var (degrees celsius)"),
    "LA_X_MIN": (0, "Lyman-alpha lower limit in x (pix)"),
    "LA_X_MAX": (0, "Lyman-alpha upper limit in x (pix)"),
    "AR_X_MIN": (0, "Active region lower limit in x (pix)"),
    "AR_X_MAX": (0, "Active region upper limit in x (pix)"),
    "AR_Y_MIN": (0, "Active region lower limit in y (pix)"),
    "AR_Y_MAX": (0, "Active region upper limit in y (pix)"),
    "PH_MIN": (0, "Minimum pulse height value (adu)"),
    "PH_MAX": (256, "Max pulse height value (adu)"),
    "PIPE_VER": (0.0, "SPRITE Pipeline Version used."),
}
