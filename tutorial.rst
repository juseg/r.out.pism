Preparing WorldClim atmospheric forcing for PISM
------------------------------------------------

This tutorial shows how to create an atmospheric forcing file for *PISM* using
the web-available *WorldClim* dataset. The resulting file contains surface
topography, temperature and precipitation data over the *European Alps*.
However, this method can be applied, with a varying degree of reality, to any
other terrestrial region north of 60°S. No prior knowledge of *GRASS GIS* is
needed. WorldClim data was choosen here because it is global and easily
processed, but it may not always be the best choice for numerical ice-sheet
modelling.

Download WorldClim data
~~~~~~~~~~~~~~~~~~~~~~~

We start by downloading altitude, temperature and precipitation data from the
WorldClim server. Go to http://www.worldclim.org/tiles.php, select tile 16 (or
whatever fits your region) and download the *altitude*, *mean temperature* and
*precipitation* files in BIL (binary array) format. You do not need to inflate
the ZIP archives. In one line for the lazy::

    wget http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/{alt,temp,prec}_16.zip

.. TODO: check file names and url

Create a new latitude/longitude location
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now let's launch GRASS GIS. Open a terminal and call::

    grass

The exact name of the program may vary depending on your version and operating
system. If you run GRASS for the first time, it will ask you to set a *GIS data
directory*. This is where all maps and temporary files will be stored. This
folder can quickly fill up a few GiB of disk space when using large datasets.

GRASS uses *locations* to organise map data by geographic projection and
coordinate system. Here we will need two locations: one corresponding to the
WorldClim data format, and another to match the desired modelling domain in
PISM.

Launch the location wizard and follow the steps to define a new location that
matches the format of WorldClim data:

* give it a comprehensive name like *world* or *world-latlon-wgs84*;
* use the *select coordinate system parameters from a list* method;
* select the coordinate system *ll* (Latitude/Longitude);
* select the ellipsoid *wgs84* (World Geodetic System 1984);
* set the default region extents to the whole globe (north=90, south=-90,
  west=-180, east=180) and leave other parameters to their default values.

Once this is done, start GRASS in the newly created location.

Import the data
~~~~~~~~~~~~~~~

Each command given below can be executed in three ways:

* by calling the command with its arguments in the terminal;
* by calling the command without arguments, which will launch a graphical
  prompt and let you enter them manually;
* by looking for the command in menus of the graphical interface and entering
  arguments manually.

To import WorldClim data into GRASS with minimal preprocessing, I invite you to
use *r.in.worldclim* (http://github.com/jsegu/r.in.worldclim/). Download the
script and run::

    r.in.worldclim.py -kyf tile=16 fields=alt,tmean,prec

Alternatively, you can inflate the downloaded files and import each layer
manually using the binary array import function *r.in.bin* from GRASS. For the
altitude field that would be::

    r.in.bin -s input=alt_16.bil output=wc_t16_alt bytes=2 north=60 south=30 \
        east=30 west=0 rows=3600 cols=3600 anull=-9999

However, WorldClim temperature and precipitation data comes in a somewhat
exotic integer format and unit conversion will be needed. It is one purpose of
*r.in.worldclim* to handle that for you.

To check if the import was successful, add the new raster layers in the layer
manager with the corresponding button in the toolbar, and refresh the display.
You should be able to visualize altitude, temperature and precipitation maps
over central Europe with a resolution of 30 arc-seconds. When you are done,
exit GRASS and restart. Since maps are stored directly on the hard disk,
nothing needs to be saved.

Create a new UTM location
~~~~~~~~~~~~~~~~~~~~~~~~~

Use the location wizard to define a new location adapted to the modelling
domain. For the Alps, we will use a UTM (Universal Transverse Mercator)
projection, zone 32, but you may need to adapt this to your domain:

* give a comprehensive name like *alps* or *alps-utm32-wgs84*;
* use the *select coordinate system parameters from a list* method;
* select the coordinate system *utm*, zone *32*;
* select the ellipsoid *wgs84* (World Geodetic System 1984);
* set the default region extents to the encompass the Alps (north=5400000,
  south=4800000, west=100000, east=1100000) and set a reasonable resolution
  (1000 m).

Start GRASS in the new location.

Reproject and export to netCDF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, reproject the topographic data in order to locate yourself::

    r.proj location=world-latlon-wgs84 input=wc_t16_alt method=bilinear

This is where GIS software becomes useful. Load the created map in the layer
manager and refresh the display. As you can see, the *default region* we have
set covers the entire Alps. This is good to produce a map, but for glacier
modelling, we can slightly reduce it in order to fit the last glacial maximum
extents. In GRASS, this is done by modifying the *current region* settings::

    g.region n=5350000 s=4850000 w=150000 e=950000

Effect a *zoom to the computational region* in the display to check if the new
bounds are reasonable. Using the same command, it is also possible to change
the resolution::

    g.region res=2000

Now that the modelling domain is defined, let's reproject the climate data. It
can be wise to reproject the topography again in order to avoid multiple
resampling::

    r.proj location=world-ll-wgs84 input=wc_t16_alt method=bilinear --overwrite
    for i in {01..12}
    do
        r.proj location=world-ll-wgs84 input=wc_t16_tmean${i} method=bilinear
        r.proj location=world-ll-wgs84 input=wc_t16_prec${i} method=bilinear
    done

Finally, export all maps with 'r.out.pism'::

    tlist=$(echo wc_t16_tmean{01..12} | tr ' ' ',')
    plist=$(echo wc_t16_prec{01..12} | tr ' ' ',')
    r.out.pism.py usurf=wc_t16_alt air_temp=${tlist} precipitation=${plist} output=atm.nc

The resulting file can be verified using for instance *ncdump* or *ncview*.
Note that *r.out.pism* automatically creates a number of variables attributes
needed by PISM, projection information, longitude and latitude maps and time
bounds. The script has several options. Please check them out using::

    r.out.pism.py --help
