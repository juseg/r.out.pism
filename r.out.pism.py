#!/usr/bin/env python

############################################################################
#
# MODULE:      r.out.pism
#
# AUTHOR(S):   Julien Seguinot
#
# PURPOSE:     Export multiple raster maps to a single NetCDF file for
#              the Parallel Ice Sheet Model [1]
#
# COPYRIGHT:   (c) 2011 Julien Seguinot
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#############################################################################

# Todo:
# * use binary tempfile when computing lon and lat
# * or generate x and y from region settings
# * improve handling of topographic maps
# * support for null values
# * nicer history line

# Note on handling of topographic maps:
# * if topg only is given:  togp=topg,                 usurf=topg;
# * if topg and thk:        topg=topg, thk=thk,        usurf=topg+thk;
# * if topg and usurf:      topg=topg, thk=usurf-topg, usurf=usurf;
# * if topg, thk and usurf: topg=topg, thk=thk,        usurf=usurf.

#%Module
#% description: Export multiple raster maps to a single NetCDF file for the Parallel Ice Sheet Model.
#% keywords: raster export netCDF PISM
#%End

#%option
#% key: topg
#% type: string
#% gisprompt: old,cell,raster
#% key_desc: name
#% description: Name of bedrock surface elevation map
#% required: yes
#%end
#%option
#% key: output
#% type: string
#% gisprompt: new,file,output
#% key_desc: name
#% description: Name for output file
#% required: yes
#%end
#%option
#% key: lon
#% type: string
#% gisprompt: old,cell,raster
#% key_desc: name
#% description: Name of longitude map (computed if lon or lat is missing)
#% required: no
#%end
#%option
#% key: lat
#% type: string
#% gisprompt: old,cell,raster
#% key_desc: name
#% description: Name of latitude map (computed if lon or lat is missing)
#% required: no
#%end
#%option
#% key: thk
#% type: string
#% gisprompt: old,cell,raster
#% key_desc: name
#% description: Name of land ice thickness map
#% required: no
#%end
#%option
#% key: usurf
#% type: string
#% gisprompt: old,cell,raster
#% key_desc: name
#% description: Name of ice surface elevation map
#% required: no
#%end
#%option
#% key: bheatflx
#% type: string
#% gisprompt: old,cell,raster
#% key_desc: name
#% description: Name of geothermal flux map
#% required: no
#%end
#%option
#% key: tillphi
#% type: string
#% gisprompt: old,cell,raster
#% key_desc: name
#% description: Name of till friction angle map
#% required: no
#%end
#%option
#% key: temp_ma
#% type: string
#% gisprompt: old,cell,raster
#% key_desc: name
#% description: Name of mean annual air temperature map
#% required: no
#%end
#%option
#% key: snowprecip
#% type: string
#% gisprompt: old,cell,raster
#% key_desc: name
#% description: Name of mean annual snow accumulation map
#% required: no
#%end
#%option
#% key: smooth
#% type: double
#% description: Neighborhood size for smoothing bed elevation
#% required: no
#%end

#%flag
#%  key: c
#%  description: Use degree Celcius instead of Kelvin
#%end

import os
import time
from numpy import *         # scientific module Numpy [2]
from netCDF4 import Dataset # interface to netCDF4 library [3]
import pyproj               # interface to PROJ.4 library [4]
from grass.script import core as grass

def readmap(mapname):	# read grass map into array
		# read current region
		region = grass.region()
		cols = int(region['cols'])
		rows = int(region['rows'])

		# read values via temporary binary file
		tempfile = grass.tempfile()
		grass.run_command('r.out.bin', input=mapname, output=tempfile, quiet=True)
		return flipud(reshape(fromfile(tempfile,dtype='f4'), (rows,cols)))
		grass.try_remove(tempfile)

def main(): # main function, called at execution time
		# parse options
		topgmap       = options['topg']       # optional
		output        = options['output']     # required
		lonmap        = options['lon']        # optional
		latmap        = options['lat']        # optional
		thkmap        = options['thk']        # optional
		usurfmap      = options['usurf']      # optional
		bheatflxmap   = options['bheatflx']   # optional
		tillphimap    = options['tillphi']    # optional
		temp_mamap    = options['temp_ma']    # optional
		snowprecipmap = options['snowprecip'] # optional
		smooth        = options['smooth']     # optional
		celcius = flags['c']

		# display warning message if both thk and usurf are defined
		if thkmap and usurfmap:
			grass.warning('thk and usurf are both defined!')

		# read current region
		region = grass.region()
		cols = int(region['cols'])
		rows = int(region['rows'])

		# ready to write NetCDF file
		ncfile = Dataset(output,'w',format='NETCDF3_CLASSIC')

		# set global attributes
		ncfile.Conventions = 'CF-1.4'
		ncfile.PROJ4Definition = grass.read_command("g.proj", flags='jf')
		ncfile.history = time.asctime() + ': ' + os.environ['CMDLINE']

		# define the dimensions
		xdim = ncfile.createDimension('x', cols)
		ydim = ncfile.createDimension('y', rows)

		# set projection coordinates
		xvar = ncfile.createVariable('x', 'f8', dimensions=('x',))
		xvar.axis = 'X'
		xvar.long_name = 'x-coordinate in Cartesian system'
		xvar.standard_name = 'projection_x_coordinate'
		xvar.units = 'm'
		for i in range(cols):
			xvar[i]=region['w']+(i-.5)*region['ewres']

		yvar = ncfile.createVariable('y', 'f8', dimensions=('y',))
		yvar.axis = 'Y'
		yvar.long_name = 'y-coordinate in Cartesian system'
		yvar.standard_name = 'projection_y_coordinate'
		yvar.units = 'm'
		for i in range(rows):
			yvar[i]=region['s']+(i-.5)*region['nsres']

		# set longitude and latitude
		lonvar = ncfile.createVariable('lon', 'f4', dimensions=('y', 'x'))
		lonvar.long_name = 'longitude'
		lonvar.standard_name = 'longitude'
		lonvar.units = 'degrees_east'

		latvar = ncfile.createVariable('lat', 'f4', dimensions=('y', 'x'))
		latvar.long_name = 'latitude'
		latvar.standard_name = 'latitude'
		latvar.units = 'degrees_north'

		if lonmap and latmap:
			# read values from grass map
			lonvar[:] = readmap(lonmap)
			latvar[:] = readmap(latmap)
		else:
			grass.message('Longitude and / or latitude map(s) unspecified. Calculating values from current projection. This can take a while...')
			# read x and y via temporary xyz file
			tempfile = grass.tempfile()
			grass.run_command('r.out.xyz', input=topgmap, output=tempfile, quiet=True)
			xy = loadtxt(tempfile, delimiter='|', dtype=int, usecols=(0,1))
			grass.try_remove(tempfile)

			# convert to longitude and latitude using pyproj
			currentproj = pyproj.Proj(grass.read_command("g.proj", flags='jf'))
			lonlat = pyproj.transform(currentproj, pyproj.Proj(proj='latlong'), xy[:,0], xy[:,1])
			lonvar[:] = flipud(reshape(lonlat[0], (rows, cols)))
			latvar[:] = flipud(reshape(lonlat[1], (rows, cols)))

		# set bedrock surface elevation
		topgvar = ncfile.createVariable('topg', 'f4', dimensions=('y', 'x'))
		topgvar.long_name = 'bedrock surface elevation'
		topgvar.standard_name = 'bedrock_altitude'
		topgvar.units = 'm'
		if smooth:
			grass.message('Smoothing bed surface elevation map...')
			newtopgmap='smoothed'
			newtopgmap = 'r.out.pism_'+str(os.getpid())+'_tmp'
			grass.run_command('r.neighbors', flags='c', input=topgmap, output=newtopgmap, size=smooth)
		grass.message('Exporting bed surface elevation map...')
		if smooth:
			topgvar[:] = readmap(newtopgmap)
			grass.run_command('g.remove', rast=newtopgmap)
		else:
			topgvar[:] = readmap(topgmap)

		# optional variable land ice thickness
		if thkmap or usurfmap:
			thkvar = ncfile.createVariable('thk', 'f4', dimensions=('y', 'x'))
			thkvar.long_name = 'land ice thickness'
			thkvar.standard_name = 'land_ice_thickness'
			thkvar.units = 'm'
			if thkvar:
				grass.message('Exporting land ice thickness map...')
				thkvar[:] = readmap(thkmap)

		# set ice surface elevation
		grass.message('Exporting ice surface elevation map...')
		usurfvar = ncfile.createVariable('usurf', 'f4', dimensions=('y', 'x'))
		usurfvar.long_name = 'ice upper surface elevation'
		usurfvar.standard_name = 'surface_altitude'
		usurfvar.units = 'm'
		if usurfmap:
			usurfvar[:] = readmap(usurfmap)
			if not thkmap:
				grass.message('Exporting land ice thickness map...')
				thkvar[:] = usurfvar[:] - topgvar[:]
		elif thkmap:
			usurfvar[:] = topgvar[:] + thkvar[:]
		else:
			usurfvar[:] = topgvar[:]

		# optional variable geothermic flux
		if bheatflxmap:
			grass.message('Exporting geothermic flux map...')
			bheatflxvar = ncfile.createVariable('bheatflx', 'f4', dimensions=('y', 'x'))
			bheatflxvar.long_name = 'upward geothermal flux at bedrock surface'
			bheatflxvar.units = 'W m-2'
			bheatflxvar[:] = readmap(bheatflxmap)

		# optional variable till friction angle
		if tillphimap:
			grass.message('Exporting till friction angle map...')
			tillphivar = ncfile.createVariable('tillphi', 'f4', dimensions=('y', 'x'))
			tillphivar.long_name = 'friction angle for till under grounded ice sheet'
			tillphivar.units = 'degrees'
			tillphivar[:] = readmap(tillphimap)

		# optional variable mean annual air temperature
		if temp_mamap:
			grass.message('Exporting mean annual air temperature map...')
			temp_mavar = ncfile.createVariable('temp_ma', 'f4', dimensions=('y', 'x'))
			temp_mavar.long_name = 'mean annual near-surface (2 m) air temperature'
			if celcius:
				temp_mavar.units = 'degC'
			else:
				temp_mavar.units = 'K'
			temp_mavar[:] = readmap(temp_mamap)

		# optional variable mean annual snow accumulation
		if snowprecipmap:
			snowprecipvar = ncfile.createVariable('snowprecip', 'f4', dimensions=('y', 'x'))
			snowprecipvar.long_name = 'mean annual ice-equivalent snow accumulation rate'
			snowprecipvar.units = 'm year-1'
			snowprecipvar[:] = readmap(snowprecipmap)

		# close NetCDF file
		ncfile.close()
		grass.message('NetCDF file '+output+' created')

if __name__ == "__main__":
    options, flags = grass.parser()
    main()

# Links
# [1] http://www.pism-docs.org
# [2] http://numpy.scipy.org
# [3] http://netcdf4-python.googlecode.com
# [4] http://pyproj.googlecode.com

