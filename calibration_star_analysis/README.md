# SPRITE CALIBRATION STAR ANALYSIS 
This set of functions inputs either one or two far or intermediate UV spectral data files and produces plots of signal-to-noise ratio per wavelength, 
photon count per pixel, photon count per wavelength, and flux per wavelength for SPRITE wavelengths at two specified exposure times.

Input:
- spectral data (.fits file)
  - from a singular mission (must be spectral data from HST, IUE, FUSE, TUES, or  HUT).
  - from two seperate missions (must be spectral data from FUSE and IUE).
  - Data should be for a standard star or some other calibration target point source.
- SPRITE effective area data (.txt file)
- SPRITE SCC effective area data (.txt file)

Output:
- 1D spectral figures showcasing the following analysis for two exposure times
  - signal-to-noise ratio per wavelength
  - photon count per pixel
  - photon count per wavelength
  - flux per wavelength

Before using this code, please make sure to change the file path in the calibration_config file.

## Summary of each .py file
"calibration_config.py" : configuration file to change variables like directory path, file path, dark rate, and more.
"get_star_data.py" :  functions to read spectral data files as well effective area data. Calls binning function.
"bin_and_convert.py" : functions to bin wavelength and flux data (2 Angstrom resolution) and convert from flux to photons.
"noise_and_calculations.py" : function to calculate signal-to-noise ratio across SPRITE wavelengths
"visualize_data.py" : functions to plot analyzed data.
"call_all_calib_funcs.py" : function to call all of the above functions.

## Demo File, Input and Output Folders

The "demo.ipynb" file is a notebook demonstrating how to run this code using input data for N132D, N49, and TwHyA as examples. (Two SNR, one Calibration Star example)

All necessary input files for the demonstration code are provided in the "sim_input" folder.

An example of simulated 2D spectral data for calibration star TwHyA is in the "sim_output" folder.


1) One data file:
  import sys
  sys.path.extend(["PATH TO YOUR CODE"])
  from call_all_calib_funcs import call_functions_one_database, call_functions_two_database
  exposure_time = [100, 1000]
  file_path = "PATH TO YOUR DATA FILES"
  file_iue = file_path  + "IUE_HD60753.fits"
  call_functions_one_database(file_iue, exposure_time, "IUE")


2) Two data files:
 
  import sys
  sys.path.extend(["PATH TO YOUR CODE"])
  from call_all_calib_funcs import call_functions_one_database, call_functions_two_database
  exposure_time = [100, 1000]
  file_path = "PATH TO YOUR DATA FILES"
  file_fuse = file_path  + "GD246_FUSE.fits"
  file_iue = file_path  + "GD246_IUE.fits"
  call_functions_two_database(file_fuse, file_iue, exposure_time)
