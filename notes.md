# Notes on lightkurve_##mag notebooks


## Notes on 7/26
I decided to export the targets table with a column with the group label for each target instead of repeating the grouping process in the notebooks.
Also, I exported the table in ECSV format which reatains the information about data type for each column.
After this changes then I'm able to simply call 
```python
targets = ascii.read('DataInput/cluster_targets_tic.ecsv')
```
to import the table with group labels, and correct data type.


Now I can remove all code in that imports and prepares the table info in each notebook. 
Also, I forgot to remove the aperture list appending after the aperture mask creation, so I got rid of it now.

Now I gotta run the corrections on each notebook.


After 14-mag the stars start becoming much less noticeable. Also, for some reason the stars look more crowded. Therefore, I decided to pump up the number of sample
TPFs that I'll use to create a master aperture mask to four TPFs. 

I noticed that about 1/3 of the stars are of 16-mag. It's probably not fair to put so many of them in a single category. 
So I'm wondering if it's worth separating them according to other factors like standard deviation, diff b/t highest value to background, etc.
I found 1/5 are in 17-mag group.

## Notes on 7/27

I made a slight modification to how I choose which pixels to keep in the master aperture masks. 
If I restart a notebook and run it again, different random files will be selected unless I seed the random generator.
Which means I accumulate more sample TPFs for the aperture mask generation. 
In order to fix that I could simply delete everything I've done, and re-run the notebooks with a **seed**.
However, as a temporary fix I opted on simply changed the way the pixels are kept for the aperture mask. 
Before, I'd keep pixels if at least two aperture masks overlaped. However, now I want to keep aperture masks if two OR 1/3 of the 
aperture masks overlap.

## Notes on 7/28

I started doing the Periodogram analysis. To setup a workflow I first need to be able to combine the quality flags with the
processed light curves. My first attempt went like this:
```python
# Open Chelsea's sample light curve, there are 5 cadences missing in downloaded raw data
quality_sample = fits.open('DataInput/ClusterQuality/Sector8_Sample.fits.gz')
# So I got to use a shorter quality mask
quality_mask = quality_sample[1].data['QUALITY'][5:]
# I replace all of the 1s with 128s, matching the quality flags: 128 = Manual Exclude
quality_mask[quality_mask == 1] = 128

# Generate TessLightCurveFile object using the saved lightcurve (FITS), and Chelsea's quality mask
# The get_lightcurve() method will output a TessLightCurve object
lc = TessLightCurveFile('DataOutput/LightCurves/TESS_LC_144752281_SEC8.fits',
                        quality_bitmask=quality_mask).get_lightcurve('FLUX')

```

However, this didn't work. The light curves that my code has generated have no quality mask. For some reason, If I try doing the above
method, the TessLightCurveFile generation will not combine my light curve's quality flags and Chelsea's `quality_mask`. Instead, I 
thought it might be a better idea to write Chelsea's quality flags into my exported light curve files, such that the quality flags
are already imported when I do:

```python
lc = TessLightCurveFile('DataOutput/LightCurves/TESS_LC_144752281_SEC8.fits').get_lightcurve('FLUX')
```

A quick solution is to go through all the files and write the quality flags ONE time. As a permanent solution, I need to implement the 
mentioned solution into the code that exports the light curves in the first place. Here is my quick solution:

```python
# Retrieve Chelsea's quality flags
quality_sample = fits.open('DataInput/ClusterQuality/Sector8_Sample.fits.gz')
quality_mask = quality_sample[1].data['QUALITY'][5:]
quality_mask[quality_mask == 1] = 128

# Gather all FITS light curve files from sector 8
fits_files = glob('DataOutput/LightCurves/TESS_LC_*_SEC8.fits')

# Loop over fits files
for fitsFile in fits_files:
    # Open fitsFile in update mode with context manager
    # The update mode allows the file to be modified on disk
    with fits.open(fitsFile, mode='update') as hdul:
        hdul[1].data['QUALITY'] = quality_mask
        # Write changes back to fitsFile
        hdul.flush()
```

After exploring the TessLightCurve & LombScarglePeriodogram I found some useful methods:

- TessLightCurve.remove_outliers : Remove outlier cadences from the light curve using sigma clipping; I like using sigma = 2.5
- TessLightCurve.to_periodogram : Use either BLS or LS to get a periodogram and find periodic behavior in light curve.
    - Some useful parameters include:
        - method='lombscargle' : use LombScargle to characterize and detect periodic signals
        - num_terms = 1 : Number of terms to use when using fourier expansion (how complex the signal is)
        - min_period = 4 : Minimum detectable period
        - max_period = 20 : Maximum detectable period
        - ov_sampling = 8
- TessLightCurve.fold : given the light curve's period, plot each period on top of each of other

## Notes on 7/29

At this point, I went proceeded to use TessLightCurve's methods to produce a clean periodogram, and plot it using TIC 144752281 as an example. 
I added vertical lines to the periodogram plot to show the detected periods. I found that **the most credible period was ~6.10 days**. *This turned out to be wrong*.
I folded the top 5 detected periodic signals to confirm that the period of 6.10 was the best. However, I later discovered that the folded light curve with the period of
6 days actually was showing 2 periods.


## Notes on 7/30

Following Stephanie's advise I plotted the top 5 credible periods for the star that I started working with (TIC: 144752281).
Then I changed the min/max period to what Stephanie suggested (0.1/28 days). Since it was pretty simple to find the peaks in the small
samples I had, increasing the sampling will need some other methodology. After some googling I found a simple function to find local maxima in signals:

```python
maxima = scipy.signal.argrelmax(mysignal)
```

Once I was able to access the maxima of the new periodogram, I plotted it, and drew red dashed vertical lines showing where the maxima values are. This is when I realized that
there were a lot of periodic signals detected of short periods < 3 days, which I think is simply noise. In order to filter those out, I decided to cut off all periods whose 
power was lower than 95% percentile. There is a total of ~15 periods after filtering. Now that the period search range is larger, I was able to see the largest signal 
was a periodic behavior of ~3 days. Once I folded the light curve using this period I noticed it looks great.


## Notes on 7/31


