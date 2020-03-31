# TESSDataExploration
Explore TESS data on open cluster IC 2391 to create custom light curves, and measure rotation periods.

Copyright 2019-2020 Jose Perez Chavez.   
This code uses Stephanie's k2spin module, and part of George Zhou download_pixfiles/getlightcurve method

Use at your own risk. This is free software made available under the MIT License. For details see the LICENSE file.

DEPENDENCIES
------------

This code makes use of the following python libraries:
    
    + Numpy
    + Scipy
    + matplotlib
    + Astropy (http://www.astropy.org/)
    + lightkurve (https://docs.lightkurve.org/)

All of them are open source and can be easily installed in any machine. 

WHAT DOES THE PIPELINE DO?
---------------------------

- All scripts perform different tasks
  - download_tpfs8: This code will use the sample TPFs included to download IC2391's 118 Full Frame cut outs. 
  Each FITS file will have the following format: TPF{ticid}_S{section}C{cutsize}.fits
  - get_TargetPSF: Create a 2x2 grid of plots showcasing the aperture selection methods for every IC2391 target.
  - getLightCurves: Extract background-subtracted light curves from all found TPFs for every aperture method, and 
  save the respective FITS files and a grid of plots with the frames and light curves per target.
  - usefulFunctions: It simply holds functions that multiple other scripts use for light curve extraction, 
  aperture selection, finding TPFs, etc.
  - getPeriodograms: This script is a work in progress and will measure the periods and output a 
  periodogram plot for every target, and make a periods vs J-K plot.