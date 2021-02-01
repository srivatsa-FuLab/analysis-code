# analysis-code
Collection of code for for data analysis 


spe_plotter.m : Plot spectra from princton-tech winspec *.spe file. 

               Additional functionality:
               - Works for both the old (blue) and new (Acton-2750) spectrometer models [Brynn_mc and NV_mc respectively]
               - Incldes manual wavelength calibration on the Acton-2750 for the 300g, 1200g and 1800g gratings (updated Jul 2020)
               - Supports automated wavelength calibration in winspec
               - Fit peaks (or dips) to lorentzian [based on user peak selection] and estimate cavity Q-factor (or ZPL linewidth)
               - Integrate area under peak (or dip)
               - save data/analysis as .png or .fig file
               
               Dependencies:
               - loadSPE.m
               
