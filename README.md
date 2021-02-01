# analysis-code
Collection of code for for data analysis 


spePlotter.m : Plot spectra from princton-tech winspec *.spe file. 

               Additional functionality:
               - Works for both old and new spectrometer models
               - Incldes manual wavelength calibration on the new spectrometer for the 300g, 1200g and 1800g gratings
               - Supports automated wavelength calibration in winspec
               - Fit peaks (or dips) to lorentzian [based on user peak selection] and estimate cavity Q-factor (or ZPL linewidth)
               - Integrate area under peak (or dip)
               - save data/analysis as .png or .fig file
               
               Dependencies:
               - loadSPE.m
               
