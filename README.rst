r.out.pism
==========

**Requires:** `GRASS GIS`_, `NetCDF4-Python`_, `Numpy`_, `PyProj`_.

**r.out.pism** is a Python script for `GRASS GIS`_ to export multiple raster maps to a NetCDF file for `PISM`_. The script exports the desired topographic, climatic and other geophysical variables together with all the necessary meta-data, and computes Cartesian and geographical coordinates using current region settings and projection parameters. This is meant to provide PISM users with an automated way of generating bootstrapping files from geographical data in any desired projection. In GRASS, type ``r.out.pism.py --help`` for a list of available options.

.. links:

.. _GRASS GIS: http://grass.osgeo.org
.. _NetCDF4-Python: http://netcdf4-python.googlecode.com
.. _NumPy: http://numpy.org
.. _PISM: http://www.pism-docs.org
.. _PyProj: http://pyproj.googlecode.com

