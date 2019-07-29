Notes on lightkurve_##mag notebooks


7/26
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

7/27

I made a slight modification to how I choose which pixels to keep in the master aperture masks. 
If I restart a notebook and run it again, different random files will be selected unless I seed the random generator.
Which means I accumulate more sample TPFs for the aperture mask generation. 
In order to fix that I could simply delete everything I've done, and re-run the notebooks with a **seed**.
However, as a temporary fix I opted on simply changed the way the pixels are kept for the aperture mask. 
Before, I'd keep pixels if at least two aperture masks overlaped. However, now I want to keep aperture masks if two OR 1/3 of the 
aperture masks overlap.

7/28

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
thought it might be a better idea to 