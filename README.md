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

All of them are open source and can be easily installed in any machine. I may be missing some.
The k2spin submodule is included in this repo; make sure to pull the submodule.

KEY FOLDERS
---------------------------
DataInput: Folder dedicated to save initial tables (e.g. IC2391 targets), and other fixed data like the cluster quality flags.

TESSCuts: This folder contains the targets Target Pixel File cutouts; downloaded using `lightkurve`.

TargetFFI: This folder showcases the comparison aperture methods used in the project for all targets in jpg's.

LightCurvesPlots: These PDFs files include a grid of plots that showcase our targets' TESS Full-Frame Images, the 
aperture methods we use, and their impact on the light curves.

DraftPeriodograms: These PDFs are simple periodogram plots produced from the three light curves for every target.


WHAT DO THE SCRIPTS DO?
---------------------------

All the details of each script's function is included in them. The following will describe the order of operation using the scripts, at least for my application. Most code will need modifications in order to work with another project specifications (e.g. another cluster).

### usefulFunctions and LCFeatureExtraction
These scripts are not really part of a order-of-operation. They simply hold useful functions that scripts use for light curve
extraction, aperture selection, finding TPFs, etc. It serves as a fixed script that can hold repetitive functions.

### download_tpfs8
This script is the first one to pay attention to. In order to get started with the photometric data you probably want to download it to your computer for easy accesss. 
The data will be available as Target Pixel Files, frame cutouts of our targets cadences.
This script will take care of it as long as you adjust the parameters indicated in it.
The script includes an alternative use that will write pre-computed apertures to the TPFs, but you will probably not use this as is.

### get_TargetPSF
This script is useful to take make a quick showcase of the aperture selection methods for every cluster (IC2391) target. 
It assumes it needs to compute the Handpicked apertures so it will need to be modified to be used; but the plotting code may be more useful.

### getLightCurves
This script is more useful than the previous for showcasing the aperture selection methods AND their impact on their respective light curves.
It will extract background-subtracted light curves from all found TPFs for every aperture method, 
and save the respective FITS files and a grid of plots with the frames and light curves per target.
This way you can easily access the light curve data later. The plots are pretty cool too.

### getPeriodograms
This script is the culmination of the previous processes. It will exract the optimal periodogram for each cluster target using the functions explored in the notebook files (discussed below). The optimal periodogram is chosen by selecting the best light curve, the one with the lowest noise. The periodogram data is saved in a FITS file including a bootstrap threshold for period validity. The periodogram is plotted and saved as well with highlighted period peaks.
The process in detail is discussed in the script, and the aperture selection is explored in the LCQualityIndicator notebook.

### masterplot
This script was used to produce preliminary period distributions, and save useful information on the processed / chosen aperture methods in a compact master table. 
This master table is used to easily access the detected periods for our cluster targets, and is used in the period_distros notebook as an example. 
While the plots produced by this script are good looking, I prefer the one from the period_distros notebook.



WHAT ARE THE NOTEBOOKS?
---------------------------

The notebooks served for testing and plotting. Not all notebooks will be useful for the reader. The key notebooks are LCQualityIndicator, PattenPeriods and period_distros. These noteoboks include plenty of documentation to follow the notebook without problem, and actually give the reader an overview of the processes of this repository.

### RemakingGetLC
This notebook started the getLightCurves script. 
It was used to optimize the background subtraction process, plotting the frame cutouts, and importing the data. 
It showcases a series of comparison between non-background subtracted and background subtracted lightcurves for both Handpicked & Percentile apertures.

### FlattenAssessment
This notebook started as a rather misconception of what the Savitzky-Golay filter is used for. I attempted to prove that I could use the filter to remove (or flatten) the trends that some light curves display. 
However, it not consistent, and in theory really shouldn't work for this type of work. 
The filter is rather too simplistic for removing the trends that our light curves show. It still has plenty of cool plots so I kept the notebook around.

### LCFeatureExtraction
A rather not-friendly notebook with no documentation. The notebook served as a starting point for creating the noise and amplitude techniques later explored. I started exploring with a Boxcar filter to estimate light curve noise, but then settled with a Gaussian filter as it best estimates the underlying signal. These concepts are implemented in LCFeatureExtraction.py.

### LCFeatureExploration
This notebook served to take a look at the distribution of amplitude-noise metrics of our light curves (three per target ~300). At the end of the day it doesn't tell me much really. It did help to see which targets suposedly have high amplitde and low noise (the second part of the notebook), but that's about it.

### LCQualityIndicator
This notebook is quite useful as it has lots of documenations on what was going through my head as I went ahead and incorporated the noise-metric to implement the automated aperture selection method. 
Since the amplitude had proved a bad indicator, not shown in notebook, I went ahead and tested the noise metric. 
The noise metric showed optimistic results for indicating low-noise light curves, and potentially useful lightcurves that the periodogram will have no problem at processing.
The notebook includes plenty of examples.

### Patten & Simon
This notebook's only purpose is to create a 'literature periods comparison'. The only other periods I found for this cluster are from Patten & Simon 1996. After running the getPeriodograms script I was able to compare the periods I got, to those from the literature and address any features / issues. I include a lot of notes on my findings in comparing these two datasets, and any complications.

### period_distros
At last this notebook served to compare our period distribution to that of modelled stellar rotators of the same age (~50 Myr) as the cluster. The modelled data comes from Matt Sean (forwarded by Stephanie Douglas). I compare two datasets with different inital conditions, but the same model. At the end of the notebook I include robust notes on unusually high periods that I found in my period distribution that cannot be trusted, and tips on improving the period detection procedure.
  

TODO
---------------------------

## Handpicked Apertures Irreproducibility
Some of the scripts try to extract a target's Handpicked aperture from the start, by running the getOrigAps function. However, this is not reproducible since it is a very manual-work and it is very specific to how I did things. Instead, we should take advantage of the Target Pixel File object's pipeline_mask property that holds an aperture that is imported from the TPF FITS file. download_tpfs8.py script writes my handpicked apertures into the downloaded TPFs, but the rest of the scripts in this repo do not use it other than getLightcurves.py. We can separately write custom apertures into our targets' TPFs, and adjust the scrits so they use those instead of running getOrigAps. 