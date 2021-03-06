# TESSDataExploration
Explore TESS data on open cluster IC 2391 to create custom light curves, and measure rotation periods.
This is a repository with all documented code produced by me, Jose Perez Chavez, as part of the research with Stephanie Douglas.
The code here showcases three aperture selection methods, their impacts on our cluster targets, and a period detection method for the TESS data.
Documentation is extensive. This README will have an overview of the script description and usage, and the notebook's purpose. 
There is a lot of in-script documentaiton as well. 
You can learn more about the project in the report written by me. 
You can view it HERE (WORK IN PROGRESS).

There is missing product files, PDFs, and images. These files are stored in my own private cloud service. 
These files include the original TESS cutouts, periodogram FITS, and master table.
I shared these product files with my advisor, [Stephanie T. Douglas](https://stephanietdouglas.science/).

Copyright 2019-2020 Jose Perez Chavez.   
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
This code uses Stephanie's k2spin (sub)module, and part of George Zhou download_pixfiles/getlightcurve method


Note on Handpicked Apertures
------------

While I am working on getting the report done, a potentially source of confusion is the handpicked apertures.
Whenever code/documentation mentions Handpicked aperture I am refering to custom apertures I designed in Summmer 2019.
After downloading the original TESS cutouts in [Target Pixel Files](https://docs.lightkurve.org/tutorials/01-target-pixel-files.html) (TPF) formats, I proceeded to determine the handpicked apertures.
The target stars were grouped by magnitude bins.
There is a total of 8 magnitude bins/groups; 10, 11, 12, 13, 14, 15, 16, 17. 
Five random targets were chosen per magnitude group, and stored in SampleTPFs folder.
An aperture that minimizes noise, and maximizes amplitde was chosen per sample target.
Each aperture is an array of the size of the TESS cutouts where a value of 3 identifies a pixel inside aperture, and 1 is outside the aperture.
If an aperture pixel is repeated twice or more among the five sample apertures, then we count it towards the final magnitude group aperture.
Finally assign the magnitude group apertures to their corresponding targets as our handpicked aperture.

The function `getOrigAps` in usefulFunctions.py performs the process mentioned above.
Unfortunately, my workflow changed later and I did not document the handpicking process. 
I permanently assigned the handpicked apertures to my the downloaded TESS cut outs (stored in products storage) as their custom aperture, `pipeline_mask`.
The property is later used in getLightCurves.py so I could access the handpicked aperture without recalculating it. 
The documentation below will give more details on any assumptions made for each script.


KEY FOLDERS
---------------------------
DataInput: Folder dedicated to save initial tables (e.g. IC2391 targets), and other fixed data like the cluster quality flags.

DataOutput/SampleTPFs: These TargetPixelFiles were handpicked to produce the handpicked apertures. The process is recorded in the usefulFuncs.getOrigAps function. There is 40 samples; 5 TPFs per G-magnitude group.

TESSCuts: This folder contains the targets Target Pixel File cutouts; downloaded using `lightkurve`.

TargetFFI: This folder showcases the comparison aperture methods used in the project for all targets in jpg's.

LightCurvesPlots: These PDFs files include a grid of plots that showcase our targets' TESS Full-Frame Images, the 
aperture methods we use, and their impact on the light curves.

DraftPeriodograms: These PDFs are simple periodogram plots produced from the three light curves for every target.

SyntheticComparison & StellarModels: This datasets/PDFs are used and exported in period_distros.ipynb. Check the notebook for details.


WHAT DO THE SCRIPTS DO?
---------------------------

All the details of each script's function is included in them. The following will describe the order of operation using the scripts, at least for my application. Most code will need modifications in order to work with another project specifications (e.g. another cluster). Also, these scrips do not accept command line parameters (yet), all parameters must be modified in the script.

### usefulFunctions and LCFeatureExtraction
These scripts are not really part of a order-of-operation. They simply hold useful functions that scripts use for light curve
extraction, aperture selection, finding TPFs, etc. It serves as a fixed module that holds repetitive functions.

### download_tpfs8
This script is your first stop at getting the data ready. 
In order to get started with the photometric data you can download it to your computer for easy accesss.
This script will download all the [Target Pixel Files](https://docs.lightkurve.org/tutorials/01-target-pixel-files.html) (TPF) of interest from your target table.
The script also contains an example of how to write a custom aperture (`pipeline_mask`) to the TPF, but will skip it since I used it for my handpicked apertures.
Once you have the appropriate paramteres simply call the script with 
```
python download_tpfs8.py
```

### get_TargetPSF
This script can make frame cutout plots of your target's TPFs with multiple apertures. 
It is useful to make quick comparisons between apertures.
It currently only works with my handpicked apertures setup, but it shows a great example on how to produce the plots in folder TargetFFI.
If handpicked aperture data is available then you can execute the script with python:
```
python get_TargetPSF.py
```


### getLightCurves
We saw the previous script can showcase different apertures for the same target.
This script will now do the same, but also show their impact on their respective processed light curves.
It will extract background-subtracted light curves from all found TPFs for every aperture method, 
and save the respective FITS files and a grid of plots with the frames and light curves per target.
This way you can easily access the light curve data later.
Now we are not assuming to have the handpicked aperture data available, but we do have custom `pipeline_mask` apertures on each TPF.
I worked on different design iterations for this plot and liked how they turned out. Feel free to use it.
Once all parameters have been tweked you can execute the script with:
```
python getLightCurves.py
```


### getPeriodograms
This script is the culmination of previous scripts and the notebooks from the next section. 
It will exract the periodogram of each target's best light curve from the three aperture methods I used; handpicked, percentile, and threshold.
The best light curve is the one that has the lowest noise.
The periodogram is saved in a FITS file including a bootstrap threshold for period validity. 
The periodogram is plotted and saved as well with highlighted potential stellar rotation periods.
The process in detail is discussed in the script, and the aperture selection is explored in the LCQualityIndicator notebook.
Once all parameters have been tweked you can execute the script with:
```
python getPeriodograms.py
```

### masterplot
This script was used to produce preliminary period distributions, and save useful information on the processed / chosen aperture methods in a compact master table. 
This master table is used to easily access the detected periods for our cluster targets, and is used in the period_distros notebook as an example. 
While the plots produced by this script are good looking, I prefer the one from the period_distros notebook. For more info check the in-script docs.
Again, for execution simply run:
```
python masterplot.py
```


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

### PattenPeriods
This notebook's only purpose is to create a 'literature periods comparison'. The only other periods I found for this cluster are from Patten & Simon 1996. After running the getPeriodograms script I was able to compare the periods I got, to those from the literature and address any features / issues. I include a lot of notes on my findings in comparing these two datasets, and any complications.

### period_distros
At last this notebook served to compare our period distribution to that of modelled stellar rotators of the same age (~50 Myr) as the cluster. The modelled data comes from Matt Sean (forwarded by Stephanie Douglas). I compare two datasets with different inital conditions, but the same model. At the end of the notebook I include robust notes on unusually high periods that I found in my period distribution that cannot be trusted, and tips on improving the period detection procedure.


TODO
---------------------------

## Handpicked Apertures Irreproducibility
The scripts included in this repository aren't very consistent when attempting to extract the handpicked aperture for each target.
download_tpfs8.py script writes my handpicked apertures into the downloaded TPFs and assigns it to pipeline_mask using getOrigAps function.
This is done as an attempt to later consistently get the handpicked apertures from the TPFs themselves. 
Furthermore, the getOrigAps function assumes my directory structure, and samples have been handpicked to then produce the apertures. 
The sample handpicking is described in the **Note on Handpicked Apertures** section above. 
The issue with this is that it assumes you already had the TPF samples per group such that we can automatically assign the handpicked apertures into each downloaded TPF's pipeline_mask property in download_tpfs8.py.
However, other scripts such as getTargetPSFs.py do not follow this rule, and attempt to extract the handpicked aperture using getOrigAps all over again.
The scripts must be modified such that there is a consistent way to get the handpicked apertures.

Altogether it becomes a bit confusing to try and reproduce the aperture handpicking. 
Which is why I stored my original TPF samples in this repo, and my downloaded TESS cutouts with the appropriate handpicked apertures are stored in my products storage. 
