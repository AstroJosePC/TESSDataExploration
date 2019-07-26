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


- 10 mag group notes
    - I re

- 11 mag group notes


- 12 mag group notes


- 13 mag group notes


- 14 mag group notes


- 15 mag group notes
    - 