# Incremental-reconstruction-of-water-tight-surface-from-unorganized-points


<h3>Requirements</h3>
VTK (version >= 6.2,https://www.vtk.org/). 
Parts of the codes are extracted from CGAL (https://www.cgal.org/).

<h3>Run</h3>

```
SurfaceReconstruction.exe data/point_cloud.xyz T
```

In release model, the command line is written like this "SurfaceReconstruction.exe data/193_r.xyz T", where 'T' denotes inserting isolated vertices, while 'F' denotes ignoring isolated vertices, and 'data/193_r.xyz' is the speicific file path. The reconstructed results will be outputted in '.stl' format, and press 's' button , a screenshot image of the model can be acquired.



