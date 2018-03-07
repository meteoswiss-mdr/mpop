"""Loader for MSG, netcdf format.
"""
from ConfigParser import ConfigParser
from mpop import CONFIG_PATH
import os
from netCDF4 import Dataset
import numpy.ma as ma
from numpy import rot90
from numpy import flipud, fliplr
from glob import glob
from mpop.projector import get_area_def
from datetime import timedelta
from numpy import logical_or

def load(satscene, **kargs):
    """Load MSG SEVIRI CPP retrieval in netCDF format.
    """

    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))
    values = {"satname": satscene.satname,
    "number": satscene.number,
    "instrument": satscene.instrument_name,
    "satellite": satscene.fullname
    }

    filename = os.path.join( satscene.time_slot.strftime(conf.get("seviri-level2", "dir", raw=True)),
                             satscene.time_slot.strftime(conf.get("seviri-level2", "filename", raw=True)) % values )

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
    ds = Dataset(filename, 'r')


    # Load specified data sets from netCDF file
    for chn_name in satscene.channels_to_load:
        # Read variable corresponding to channel name
        data = ds.variables[chn_name][:,:]   #
        data = rot90(data, k=3)
        #data = flipud(data)
        #data = fliplr(data)
        print '(A min/max) ',data.min(),data.max()

        #cot:units = "1" ;
        try:
            units = ds.variables[chn_name].getncattr("units")
            if units=="1":
                units=None
            print "... units ", units
        except AttributeError:
            units=None

        try:
            FillValue = ds.variables[chn_name].getncattr("_FillValue")
            print "... FillValue ", FillValue
            data = ma.masked_equal( data, FillValue )
            data = ma.masked_equal( data, 0 )
        except AttributeError:
            FillValue=None

        try:
            long_name = ds.variables[chn_name].getncattr("long_name")
            print "... long_name ", long_name
        except AttributeError:
            long_name=None

        print type(data)

        #from trollimage.image import Image as trollimage
        #img = trollimage(data[::3,::3], mode="L") # [0,0,0] [1,1,1], fill_value=[1,1,1]
        #from trollimage.colormap import rainbow
        #img.colorize(rainbow)
        #img.show()
        #quit()


        print '(B min/max) ',data.min(),data.max()

        satscene[chn_name] = ma.asarray(data)
        print "*** *** units ", units
        if units != None:
            satscene[chn_name].info['units'] = units

        # !!! BAD: THIS INFORMATION SHOULD BE SAVED IN THE FILE !!!
        satscene.area = get_area_def("SeviriDiskFull00")

        #satscene.channels.append(ct_chan)
