# SPRITE_DRP
SPRITE Data Reduction Pipeline and Sample Data Products

* **mast_sample_files/** contains the scripts to generate sample SPRITE data products, which are then stored in the 'generated/' subdirectory. These files have no actual data - they are just the data structures and metadata in the headers. 
* **simulated_observations/** contains Elena's code and data for simulating SPRITE observations images. These don't have proper headers or anything - just the simulated images.
* **drp_main/** contains the data reduction routines and standard star data to be used for flux & wavelength calibration.