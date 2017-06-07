"""Loader for MSG, netcdf format.
"""
from ConfigParser import ConfigParser
from mpop import CONFIG_PATH
import os
##from netCDF4 import Dataset  ### !!! dont use this !!!
from Scientific.IO.NetCDF import NetCDFFile
import numpy.ma as ma
from glob import glob
from mpop.projector import get_area_def
from numpy import flipud, ndarray

long_names = {'vza': 'viewing zenith angle', 'vaa': 'viewing azimuth angle', 'lat': 'latitude',     'lon': 'longitude'}
units      = {'vza': 'degree',               'vaa': 'degree',                'lat': 'degree north', 'lon': 'degree east'}

def load(satscene, **kargs):
    """Load MSG SEVIRI viewing geometry from netCDF file.
    """

    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))
    values = {"orbit": satscene.orbit,
    "satname": satscene.satname,
    "number": satscene.number,
    "instrument": satscene.instrument_name,
    "satellite": satscene.fullname
    }


    # get substring between 'lon_0=' and the following space  
    proj4_params = conf.get("satellite", "proj4_params", raw=True)
    lon_0 = float(proj4_params [ proj4_params.find('lon_0')+6: proj4_params.find(' ',proj4_params.find('lon_0')) ])
    # print "lon_0 ", lon_0
    lon0_str = '{0:02d}'.format(int(10*lon_0))
    values["lon_0"] = lon0_str

    filename = os.path.join( satscene.time_slot.strftime(conf.get("seviri-level6", "dir", raw=True)),
                             satscene.time_slot.strftime(conf.get("seviri-level6", "filename", raw=True)) % values )

    print "... search for file: ", filename
    filenames=glob(str(filename))
    if len(filenames) == 0:
        print "*** Error, no file found"
        quit()
    elif len(filenames) > 1:
        print "*** Warning, more than 1 datafile found: ", filenames 
    filename = filenames[0]
    print("... read data from %s" % str(filename))

    # Load data from netCDF file
    ncfile = NetCDFFile(filename, 'r')
    for chn_name in satscene.channels_to_load: 
        # Read variable corresponding to channel name
        data = ncfile.variables[chn_name][:]   # attention [:,:] or [:] is really necessary
        print type(data)                       # it converts 'NetCDFVariable' to 'numpy.ndarray'
        print data.ndim
        print data.shape

        data_m = ma.masked_array(data) # convert numpy array to masked array 
        data_m.mask = (data == 999.0)  # create mask 

        # with this rotation, we can use same indices as we do in the EUMETSAT order mask 
        # data_m.mask[3221:3371, 1712:2060] = 1

        #import matplotlib.pyplot as plt
        #plt.imshow(data_m)
        #plt.colorbar()
        #plt.show()

        satscene[chn_name] = flipud(data_m)   # flip upside-down 
        satscene[chn_name].fill_value = 999.0
        satscene[chn_name].info['units'] = units[chn_name]
        satscene[chn_name].long_name  = long_names[chn_name]

        # extract string for the area from filename 
        area = filename[filename.rfind("_")+1:filename.rfind(".")]
        print "... set area to ", area 
        # !!! BAD: THIS INFORMATION SHOULD BE SAVED IN THE FILE !!!
        satscene.area = get_area_def(area)

    # close the file.
    ncfile.close()
