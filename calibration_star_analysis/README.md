This set of functions inputs either one or two FUV spectral data files and produces plots of signal-to-noise ratio per wavelength, 
photon count per pixel, photon count per wavelength, and flux per wavelength.

Inputs include single files (must be spectral data from HST, IUE, FUSE, TUES, or  HUT) as well as two files (must be spectral data from FUSE and IUE).  Data should
be for a standard star or some other calibration target point source.  Other input files include effective area curve data, SCC effective 
area curve data.

Before using this code, please make sure to change the file path in the calibration_config file.

Below are two ways to call the functions:

1) One data file:
  import sys
  sys.path.extend(["*** PATH TO YOUR CODE ****"])
  from call_all_calib_funcs import call_functions_one_database, call_functions_two_database
  exposure_time = [100, 1000]
  file_path = "*** PATH TO YOUR DATA FILES ***"
  file_iue = file_path  + "IUE_HD60753.fits"
  call_functions_one_database(file_iue, exposure_time, "IUE")


2) Two data files:
 
  import sys
  sys.path.extend(["*** PATH TO YOUR CODE ****"])
  from call_all_calib_funcs import call_functions_one_database, call_functions_two_database
  exposure_time = [100, 1000]
  file_path = "*** PATH TO YOUR DATA FILES ***"
  file_fuse = file_path  + "GD246_FUSE.fits"
  file_iue = file_path  + "GD246_IUE.fits"
  call_functions_two_database(file_fuse, file_iue, exposure_time)
