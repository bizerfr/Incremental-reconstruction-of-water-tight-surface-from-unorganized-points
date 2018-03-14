# Incremental-reconstruction-of-water-tight-surface-from-unorganized-points

This code is used to reconstruct water-tight and genus-0 surface from points. 
See "http://intlpress.com/site/pub/pages/journals/items/cis/content/vols/0017/0002/a002/index.html" for speicific details.
To implement it, VTK is required. 
Some parts of the code are refereced by CGAL.
For better performance, rearrange the order of '.xyz' file randomly. You could use the 'xyz_random.m' to achieve this goal.

In release model, the command line is written like this "SurfaceReconstruction.exe data/193_r.xyz T", where 'T' denotes inserting isolated vertices, while 'F' denotes ignoring isolated vertices, and 'data/193_r.xyz' is the speicific file path.



