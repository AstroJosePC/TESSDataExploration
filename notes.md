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

I made a slight modification to how I choose which pixels to keep in the master aperture masks. 
If I restart a notebook and run it again, different random files will be selected unless I seed the random generator.
Which means I accumulate more sample TPFs for the aperture mask generation. 
In order to fix that I could simply delete everything I've done, and re-run the notebooks with a **seed**.
However, as a temporary fix I opted on simply changed the way the pixels are kept for the aperture mask. 
Before, I'd keep pixels if at least two aperture masks overlaped. However, now I want to keep aperture masks if two OR 1/3 of the 
aperture masks overlap.

- 10 mag group notes
    - I re

- 11 mag group notes


- 12 mag group notes


- 13 mag group notes


- 14 mag group notes


- 15 mag group notes
    - 