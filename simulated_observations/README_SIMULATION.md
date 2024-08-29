### SPRITE SIMULATED OBSERVATIONS

This code creates 1D and 2D spectral simulations that predict what the SPRITE CubeSat will observe when looking at different supernova remnants (SNR) or calibration stars in the Large Magellanic Cloud (LMC) and Small Magallenic Cloud (SMC).  Produces a low resolution image of the reconstructed SNR.

Input files:
- SPRITE detector data (.fits file)
- FUSE and IUE SNR spectra (.fits files)
- SPRITE effective area data (.txt file)

Output:
- 1D spectral model (simulated)
- 2D spectral model with noise (simulated)


## Summary of each .py file:

1. "load_supernova_data.py" - functions to read in FUSE and IUE .fits files as well as SPRITE detector data.  Calls binning function and combines FUSE and IUE data to span SPRITE wavebands.  Calls conversion function to convert flux values to photon counts.

2. "bin_and_convert_data.py" - functions to bin wavelength and flux data (1 Angstrom resolution) and convert from flux to photons

3. "create_2d_supernova_spectra.py" - functions to create a 2D array of SNR data convolved onto SPRITE slit

4. "supernova_spec_noise.py" - functions to incorportae shot noise, dark current background, and geocoronal lyman alpha emission noise into the 2D spectrum

5. "reconstruct_supernova_image.py" - functions to combine brightness profiles along pushbroom steps to reconstruct a low resolution image of the SNR.

6. "plot_supernova_data.py" - functions to plot 1D and 2D simulated spectra, as well as reconstructed SNR image

7. "sprite_sim_config.py" - configuration file to change variables like directory paths, exposure time, detector size, dark rate, and more.


## Demo File, Input and Output Folders
The "demo.ipynb" file is a notebook demonstrating how to run this code using N132D input data as an example.

All necessary input files for the demonstration code are provided in the "sim_input" folder.

An example of simulated 2D spectral data for calibration star TwHyA is in the "sim_output" folder.





