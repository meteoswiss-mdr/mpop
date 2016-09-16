"""Loader for MSG, netcdf format.
"""
from ConfigParser import ConfigParser
from mpop import CONFIG_PATH
import os
from netCDF4 import Dataset
import numpy.ma as ma
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

    time_slot2 = satscene.time_slot + timedelta(minutes=15)
    enddate = time_slot2.strftime("%Y%m%dT%H%M00")
    print enddate
    values["enddate"] = enddate

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

    # read solar zenith angle to create data mask 
    sunz = ds.variables['sunz'][0,:,:]       # !!! change this !!! take first time slice 
    satz = ds.variables['satz'][0,:,:]       # !!! change this !!! take first time slice 
    mask = ma.masked_array( logical_or(sunz > 78.0, satz > 78.0) )
    mask = ma.masked_equal( mask, False)
    satscene['MASK'] = mask

    #from trollimage.image import Image as trollimage
    #img = trollimage(mask, mode="L", fill_value=[1,1,1]) # [0,0,0] [1,1,1]
    #from trollimage.colormap import rainbow
    #img.colorize(rainbow)
    #img.show()
    #quit()

    # Load specified data sets from netCDF file
    for chn_name in satscene.channels_to_load:
        # Read variable corresponding to channel name
        data = ds.variables[chn_name][0,:,:]   # !!! change this !!! take first time slice 

        #print '(A min/max) ',data.min(),data.max()

        #cot:units = "1" ;
        units = ds.variables[chn_name].getncattr("units")
        if units=="1":
            units=None
        print "... units ", units

        if chn_name == 'cph':
            data = ma.masked_equal( data, 0 ) # 0 == no clouds
        if chn_name == 'cot':
            data = ma.masked_equal( data, 0 ) # 0 == no clouds
        if chn_name == 'cth':
            data = data / 1000. # m -> km
            units = 'km'
        if chn_name == 'dcld':
            data = data / 1000. # m -> km
        #if chn_name == 'dndv':
        #    data = data / 1.E9  # particle / m3 -> particles / mm3
        if chn_name == 'dreff':
            data = data * 1.E6  # m -> micro meter
            units = 'micro m'
        if chn_name == 'precip' or chn_name == 'precip_ir':
            data = (1000 * 3600) * data   # m / s -> 1000mm / 3600s 
            data = ma.masked_equal( data, 0 ) # 0 == no precip
            units = 'mm/h'
        if chn_name == 'reff':
            data = data * 1.E6  # m -> micro meter
            data = ma.masked_equal( data, 0 ) # 0 == no clouds
            units = 'micro m'

        #print '(B min/max) ',data.min(),data.max()

        satscene[chn_name] = ma.asarray(data)
        print "*** *** units ", units
        if units != None:
            satscene[chn_name].info['units'] = units
        # !!! BAD: THIS INFORMATION SHOULD BE SAVED IN THE FILE !!!
        satscene.area = get_area_def("SeviriDiskFull00")

        #satscene.channels.append(ct_chan)
