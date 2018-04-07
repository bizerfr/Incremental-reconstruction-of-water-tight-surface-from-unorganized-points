# Incremental-reconstruction-of-water-tight-surface-from-unorganized-points

This code is used to reconstruct water-tight and genus-0 surface from points. 
See "http://intlpress.com/site/pub/pages/journals/items/cis/content/vols/0017/0002/a002/index.html" and "http://kns.cnki.net/KCMS/detail/51.1196.TP.20171010.1731.106.html?uid=WEEvREcwSlJHSldRa1Fhb09jSnZqYVRDOElYcUw4ZERvYzAxU29TUVp0VT0=$9A4hF_YAuvQ5obgVAqNKPCYcEjKensW4ggI8Fm4gTkoUKaID8j8gFw!!&v=MDI0MjZyTUpGWT1MejdTWkxHNEg5Yk5yNDVGWk9weVl3OU16bVJuNmo1N1QzZmxxV00wQ0xMN1I3cWViK1pyRmk3bVU3" for speicific details.

To implement it, VTK (version >= 6.2,https://www.vtk.org/) is required. 
Some algorithmic parts of the code are refereced by CGAL(https://www.cgal.org/, implementation of this code is free installation of CGAL).
For better performance, rearrange the order of '.xyz' file randomly. You could use the 'xyz_random.m' to achieve this goal.

In release model, the command line is written like this "SurfaceReconstruction.exe data/193_r.xyz T", where 'T' denotes inserting isolated vertices, while 'F' denotes ignoring isolated vertices, and 'data/193_r.xyz' is the speicific file path. The reconstructed results will be outputted in '.stl' format, and press 's' button , a screenshot image of the model can be acquired.



