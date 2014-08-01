#!/usr/bin/env python2

"""
MODULE:     r.out.pism

AUTHOR(S):  Julien Seguinot

PURPOSE:    Export multiple raster maps to a single NetCDF file for
            the Parallel Ice Sheet Model [1]

COPYRIGHT:  (c) 2011-2014 Julien Seguinot

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# Todo:
# * add support for null values
# * add option for arbitrary order of dimensions
# * use grass.array

# Version history
# * 14/01/2013 (v0.3)
#  - added bytes=4 in export command to prevent datatypes error
# * 19/09/2012
#  - added time_bounds variable
# * 05/07/2012
#  - added edgetemp option
#  - corrected erroneous offset of x and y variables by one grid space
# * 01/06/2012
#  - updated variable names for pism0.5 compatibility
# * 28/09/2011 (v0.2r1)
#  - added artm for compatibility with bc files
#  - added mapping variable with proj.4 attribute for use in PISM
# * 01/07/2011 (v0.2)
#  - changed order of netCDF dimensions to txyz
#  - added time dimension and multiple temp_ma and precip export
# * 13/06/2011:
#  - changed snowprecip to precip
# * 27/05/2011:
#  - added -f flag to use Fahrenheit degrees
#  - added -x flag to turn off lon / lat calculation
#  - applied optional smoothing to all maps instead of just topg
#  - fastened the lon / lat calculation (no tempfile)
#  - customized the handling of topographic maps
# * 13/05/2011 (v0.1):
#  - first version of the script

# Note on the handling of topographic maps:
# |input \ output      |topg      |thk        |usurf
# |====================================================
# |topg only           |topg      | -         | -
# |thk only            | -        |thk        | -
# |usurf only          | -        | -         |usurf
# |topg and thk        |topg      |thk        |topg+thk
# |topg and usurf      |topg      |usurf-topg |usurf
# |thk and usurf       |usurf-thk |thk        |usurf
# |topg, thk and usurf |topg      |thk        |usurf

#%Module
#% description: Export multiple raster maps to a NetCDF file for PISM.
#% keywords: raster export netCDF PISM
#%End

#%option
#% key: output
#% type: string
#% gisprompt: new,file,output
#% description: Name for output file
#% required: yes
#%end
#%option
#% key: lon
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of longitude map (computed if lon or lat is missing)
#% required: no
#%end
#%option
#% key: lat
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of latitude map (computed if lon or lat is missing)
#% required: no
#%end
#%option
#% key: topg
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of bedrock surface elevation map
#% required: no
#% multiple: yes
#%end
#%option
#% key: thk
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of land ice thickness map
#% required: no
#% multiple: yes
#%end
#%option
#% key: usurf
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of ice surface elevation map
#% required: no
#% multiple: yes
#%end
#%option
#% key: bheatflx
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of geothermal flux map
#% required: no
#% multiple: yes
#%end
#%option
#% key: tillphi
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of till friction angle map
#% required: no
#% multiple: yes
#%end
#%option
#% key: air_temp
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of near-surface air temperature map
#% required: no
#% multiple: yes
#%end
#%option
#% key: air_temp_sd
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of standard deviation of near-surface air temperature map
#% required: no
#% multiple: yes
#%end
#%option
#% key: temp_ma
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of mean annual air temperature map
#% required: no
#% multiple: yes
#%end
#%option
#% key: precipitation
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of mean annual snow accumulation map
#% required: no
#% multiple: yes
#%end
#%option
#% key: edgetopg
#% type: double
#% description: Bedrock surface elevation at the edge of the domain
#% required: no
#% answer: -3000
#%end
#%option
#% key: edgetemp
#% type: double
#% description: Temperature value for the edge of the domain
#% required: no
#% answer: 300
#%end
#%option
#% key: smooth
#% type: integer
#% description: Neighborhood size for optional smoothing before import
#% required: no
#%end

#%flag
#%  key: c
#%  description: Use degree Celcius instead of Kelvin
#%end
#%flag
#%  key: e
#%  description: Add buffer zone along edges to prevent ice leak.
#%end
#%flag
#%  key: f
#%  description: Use degree Fahrenheit instead of Kelvin
#%end
#%flag
#%  key: p
#%  description: Divide precipitation by ice density (0.910)
#%end
#%flag
#%  key: x
#%  description: Do not compute lon and lat if missing.
#%end

import os
import sys
import time
from numpy import *                     # scientific module Numpy [2]
from netCDF4 import Dataset, Variable   # interface to netCDF4 library [3]
import pyproj                           # interface to PROJ.4 library [4]
from grass.script import core as grass
from grass.script import array as garray


### Default variable names and attributes ###

names = {

    # coordinate variables
    'x': {
        'axis': 'X',
        'long_name': 'x-coordinate in Cartesian system',
        'standard_name': 'projection_x_coordinate',
        'units': 'm'},
    'y': {
        'axis': 'Y',
        'long_name': 'y-coordinate in Cartesian system',
        'standard_name': 'projection_y_coordinate',
        'units': 'm'},
    'lon': {
        'long_name': 'longitude',
        'standard_name':'longitude',
        'units': 'degrees_east'},
    'lat': {
        'long_name': 'latitude',
        'standard_name':'latitude',
        'units': 'degrees_north'},
    'time': {
        'axis': 'T',
        'long_name': 'time',
        'standard_name': 'time',
        'bounds': 'time_bounds',
        'units': 'month'},
    'time_bounds': {},
    'mapping': {},

    # topographic variables
    'topg': {
        'long_name': 'bedrock surface elevation',
        'standard_name': 'bedrock_altitude',
        'units': 'm'},
    'thk': {
        'long_name': 'land ice thickness',
        'standard_name': 'land_ice_thickness',
        'units': 'm'},
    'usurf': {
        'long_name': 'ice upper surface elevation',
        'standard_name': 'surface_altitude',
        'units': 'm'},

    # climatic variables
    'air_temp': {
        'long_name': 'near-surface air temperature',
        'units': 'K'},
    'air_temp_sd': {
        'long_name': 'standard deviation of near-surface air temperature',
        'units': 'K'},
    'precipitation': {
        'long_name': 'mean annual %sprecipitation rate',
        'standard_name': '',
        'units': 'm year-1'},
    'temp_ma': {
        'long_name': 'mean annual near-surface (2 m) air temperature',
        'units': 'K'},

    # other variables
    'bheatflux': {
        'long_name': 'upward geothermal flux at bedrock surface',
        'units': 'mW m-2'},
    'tillphi': {
        'long_name': 'friction angle for till under grounded ice sheet',
        'units': 'degrees'}}


### Internal functions ###

def grass_str_list(option):
    """Return a list of strings from grass parsed option"""
    return option.split(',') if option else []


def get_dim(maplist):
    """Return dimension tuple to be used in NetCDF file"""
    if len(maplist) == 1:
        return ('x', 'y')
    else:
        return ('time', 'x', 'y')


def read_map(mapname, scalefactor=1.0):
    """Return numpy array from a GRASS raster map"""

    # show which map is processed if verbose
    grass.verbose(mapname)

    # parse smoothing option
    smooth = options['smooth']

    # smooth map with r.neighbors
    if smooth:
        smoothmap = 'r.out.pism_'+str(os.getpid())+'_tmp'
        grass.run_command('r.neighbors', flags='c',
                          input=mapname, output=smoothmap,
                          size=options['smooth'], quiet=True)
        mapname = smoothmap

    # read map into array
    a = garray.array()
    a.read(mapname)
    if smooth:
        grass.run_command('g.remove', rast=smoothmap, quiet=True)
    return transpose(flipud(a[:]))*scalefactor


### Customized NetCDF classes ###

class PISMDataset(Dataset):
    def createVariable(self, varname, datatype, dimensions=()):
        """Create a new variable and set default attributes"""
        var = PISMVariable(self, varname, datatype, dimensions)
        for (attr, value) in names[varname].iteritems():
            setattr(var, attr, value)
        self.variables[varname] = var
        return self.variables[varname]


class PISMVariable(Variable):
    def set_maps(self, maplist, scalefactor=1.0):
        """Set a list of GRASS raster maps as variable data"""

        # count number of maps to import
        nmaps = len(maplist)

        # if only one map, import it as 2D data
        if nmaps == 1:
            self[:] = read_map(maplist[0], scalefactor=scalefactor)

        # if several maps, import them as time slices
        else:
            for i in range(nmaps):
                self[i] = read_map(maplist[i], scalefactor=scalefactor)


### Main function ###

def main():
    """main function, called at execution time"""
    # parse arguments
    output = options['output']
    lonmap = options['lon']
    latmap = options['lat']
    topg = grass_str_list(options['topg'])
    thk = grass_str_list(options['thk'])
    usurf = grass_str_list(options['usurf'])
    bheatflx = grass_str_list(options['bheatflx'])
    tillphi = grass_str_list(options['tillphi'])
    air_temp = grass_str_list(options['air_temp'])
    air_temp_sd = grass_str_list(options['air_temp_sd'])
    temp_ma = grass_str_list(options['temp_ma'])
    precipitation = grass_str_list(options['precipitation'])
    edgetemp = options['edgetemp']
    edgetopg = options['edgetopg']
    iceprecip = flags['p']
    celcius = flags['c']
    edgebuffer = flags['e']
    fahrenheit = flags['f']
    nolonlat = flags['x']

    # multiple input is not implemented for topographic maps
    if len(topg) > 1:
        grass.fatal('Multiple topg export is not implemented yet, sorry.')
    if len(thk) > 1:
        grass.fatal('Multiple thk export is not implemented yet, sorry.')
    if len(usurf) > 1:
        grass.fatal('Multiple usurf export is not implemented yet, sorry.')

    # this is here until order of dimensions becomes an option
    twodims = ('x', 'y')

    # read current region
    region = grass.region()
    cols = int(region['cols'])
    rows = int(region['rows'])

    # read current projection
    proj = pyproj.Proj(grass.read_command("g.proj", flags='jf'))

    # open NetCDF file
    nc = PISMDataset(output, 'w', format='NETCDF3_CLASSIC')

    # set global attributes and projection info
    nc.Conventions = 'CF-1.4'
    nc.history = '%s: %s' % (time.asctime(), ' '.join(sys.argv))
    nc.proj4 = proj.srs.rstrip()
    mapping = nc.createVariable('mapping', byte)
    mapping.proj4 = proj.srs.rstrip()

    # define the dimensions
    nc.createDimension('time', None)  # None means unlimited
    nc.createDimension('x', cols)
    nc.createDimension('y', rows)
    nc.createDimension('nv', 2)

    # set projection x coordinate
    xvar = nc.createVariable('x', 'f8', ('x',))
    for i in range(cols):
        xvar[i] = region['w'] + (i+.5)*region['ewres']

    # set projection y coordinate
    yvar = nc.createVariable('y', 'f8', ('y',))
    for i in range(rows):
        yvar[i] = region['s'] + (i+.5)*region['nsres']

    # initialize longitude and latitude
    if (lonmap and latmap) or not nolonlat:
        lonvar = nc.createVariable('lon', 'f4', twodims)
        latvar = nc.createVariable('lat', 'f4', twodims)

    # export lon and lat maps if both available
    if lonmap and latmap:
        grass.message('Exporting longitude...')
        lonvar[:] = read_map(lonmap)
        grass.message('Exporting latitude...')
        latvar[:] = read_map(latmap)

    # else compute them using pyproj
    elif not nolonlat:
        grass.message('Longitude and / or latitude map(s) unspecified.'
                      'Calculating values from current projection...')
        x = repeat(xvar, rows)
        y = tile(yvar, cols)
        lonvar[:], latvar[:] = proj(x, y, inverse=True)

    # initialize bedrock surface elevation
    if topg or (thk and usurf):
        topgvar = nc.createVariable('topg', 'f4', twodims)

    # initialize land ice thickness
    if thk or (topg and usurf):
        thkvar = nc.createVariable('thk', 'f4', twodims)

    # initialize ice surface elevation
    if usurf or (topg and thk):
        usurfvar = nc.createVariable('usurf', 'f4', twodims)

    # export available topographic maps
    if topg:
        grass.message('Exporting bed surface elevation...')
        topgvar.set_maps(topg)
    if thk:
        grass.message('Exporting land ice thickness...')
        thkvar.set_maps(thk)
    if usurf:
        grass.message('Exporting ice surface elevation...')
        usurfvar.set_maps(usurf)

    # possibly compute the rest
    if not topg and (thk and usurf):
        grass.message('Computing land ice thickness...')
        topgvar[:] = usurfvar[:] - thkvar[:]
    if not thk and (topg and usurf):
        grass.message('Computing land ice thickness...')
        thkvar[:] = usurfvar[:] - topgvar[:]
    if not usurf and (topg and thk):
        grass.message('Computing ice surface elevation...')
        usurfvar[:] = topgvar[:] + thkvar[:]

    # assign given edge topography at domain edges
    if (topg or (thk and usurf)) and edgebuffer:
        topgvar[:3, :] = topgvar[-3:, :] = edgetopg
        topgvar[:, :3] = topgvar[:, -3:] = edgetopg

    # set geothermic flux
    if bheatflx:
        bheatflxvar = nc.createVariable('bheatflx', 'f4', twodims)
        grass.message('Exporting geothermic flux...')
        bheatflxvar.set_maps(bheatflx)

    # set till friction angle
    if tillphi:
        tillphivar = nc.createVariable('tillphi', 'f4', twodims)
        grass.message('Exporting till friction angle...')
        tillphivar.set_maps(tillphi)

    # set near-surface air temperature (air_temp)
    if air_temp:
        air_tempvar = nc.createVariable('air_temp', 'f4', get_dim(air_temp))
        if celcius:
            air_tempvar.units = 'degC'
        elif fahrenheit:
            air_tempvar.units = 'degF'
        grass.message('Exporting near-surface air temperature...')
        air_tempvar.set_maps(air_temp)

    # assign given edge temperature at domain edges
    if air_temp and edgebuffer:
        for i in range(air_tempvar.shape[0]):
            air_tempvar[i, :3, :] = air_tempvar[i, -3:, :] = edgetemp
            air_tempvar[i, :, :3] = air_tempvar[i, :, -3:] = edgetemp

    # set standard deviation of near-surface air temperature (air_temp_sd)
    if air_temp_sd:
        air_temp_sdvar = nc.createVariable(
            'air_temp_sd', 'f4', get_dim(air_temp_sd))
        grass.message('Exporting standard deviation of near-surface air'
                      'temperature...')
        air_temp_sdvar.set_maps(air_temp_sd)

    # set mean annual air temperature (temp_ma)
    if temp_ma:
        temp_mavar = nc.createVariable('temp_ma', 'f4', get_dim(temp_ma))
        if celcius:
            temp_mavar.units = 'degC'
        elif fahrenheit:
            temp_mavar.units = 'degF'
        grass.message('Exporting mean annual air temperature...')
        temp_mavar.set_maps(temp_ma)

    # set annual snow precipitation
    if precipitation:
        precipitationvar = nc.createVariable(
            'precipitation', 'f4', get_dim(precipitation))
        precipitationvar.long_name=('mean annual %sprecipitation rate'
                                 % ('ice-equivalent ' if iceprecip else ''))
        grass.message('Exporting precipitation rate...')
        precipitationvar.set_maps(
            precipitation,
            scalefactor=(1/0.91 if iceprecip else 1.0))

    # set time coordinate and time bounds
    timevar = nc.createVariable('time', 'f8', ('time',))
    time_boundsvar = nc.createVariable('time_bounds', 'f8', ('time', 'nv'))
    for i in range(len(nc.dimensions['time'])):
        timevar[i] = i
        time_boundsvar[i, :] = [i, i+1]

    # close NetCDF file
    nc.close()
    grass.message('NetCDF file '+output+' created')

### Main program ###

if __name__ == "__main__":
    options, flags = grass.parser()
    main()

# Links
# [1] http://www.pism-docs.org
# [2] http://numpy.scipy.org
# [3] https://github.com/Unidata/netcdf4-python
# [4] http://pyproj.googlecode.com
# [5] http://www.pism-docs.org/doxy/html/std_names.html
