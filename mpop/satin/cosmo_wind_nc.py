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
from numpy import logical_or, flipud, array, stack, append

def get_area_for_comsmo(filename):
    """Get the area from the h5 file.
    """
    from pyresample.geometry import AreaDefinition

    aex = get_area_extent(filename)
    h5f = h5py.File(filename, 'r')
    pname = h5f.attrs["PROJECTION_NAME"]
    proj = {}
    if pname.startswith("Swiss coordinate system"):
        proj["proj"] = "somerc"
        proj["ellps"] = "bessel"
        proj["k_0"] = "1"
        proj["lat_0"] = "46.9524055555556"
        proj["lon_0"] = "7.43958333333333"
        #proj["lon_0"] = str(float(pname.split("<")[1][:-1]))
    else:
        raise NotImplementedError("Only geos projection supported yet.")

    #h5f.attrs["REGION_NAME"]  # <type 'numpy.string_'> alps
    #pname                     # <type 'numpy.string_'> GEOS<+009.5>
    #proj                      # <type 'dict'> {'a': '6378169.0', 'h': '35785831.0', 'b': '6356583.8', 'lon_0': '9.5', 'proj': 'geos'}               
    #int(h5f.attrs["NC"])      # <type 'int'>  349
    #int(h5f.attrs["NL"])      # <type 'int'> 151
    #aex                       # <type 'tuple'> (-613578.17189778585, 4094060.208733994, 433553.97518292483, 4547101.2335793395)

    area_def = AreaDefinition(h5f.attrs["REGION_NAME"],
                              h5f.attrs["REGION_NAME"],
                              pname,
                              proj,
                              int(h5f.attrs["NC"]),
                              int(h5f.attrs["NL"]),
                              aex)
    h5f.close()
    return area_def


def load(satscene, **kwargs):
    """Load COSMO analysis and forecasts from netCDF format.
    """

    if "forecast_time" in kwargs:
        key = "forecast_time"
        forecast_time = kwargs[key]
        print "... %s: %s" % (key, kwargs[key])
    else:
        forecast_time = [0]

    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))
    values = {"satname": satscene.satname,
    "number": satscene.number,
    "instrument": satscene.instrument_name,
    }

    for ftime in forecast_time:

        ftime_str = "%02d" % ftime # in hours
        values["forecast_time"] = ftime_str
        print '... forecast time: '+ftime_str

        filename = os.path.join( satscene.time_slot.strftime(conf.get("cosmo-level2", "dir", raw=True)),
                                 satscene.time_slot.strftime(conf.get("cosmo-level2", "filename", raw=True)) % values )

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

        # read the dimensions
        x_1 = ds.variables['x_1']
        y_1 = ds.variables['y_1']
        z_1 = ds.variables['z_1']
        satscene.x_1 = x_1[:]
        satscene.z_1 = y_1[:]
        satscene.z_1 = z_1[:]
        satscene.x_1_long_name = ds.variables['x_1'].getncattr("long_name")
        satscene.y_1_long_name = ds.variables['y_1'].getncattr("long_name")
        satscene.z_1_long_name = ds.variables['z_1'].getncattr("long_name")

        #lon_1 = ds.variables['lon_1']
        #lat_1 = ds.variables['lat_1']

        # !!! BAD: THIS INFORMATION SHOULD BE SAVED IN THE FILE !!!
        satscene.area = get_area_def("ccs4")

        if satscene.x_1_long_name != "x-coordinate in Swiss coordinate system" :
            print "*** Error in mpop/satin/cosmo2_wind_nc.py"
            print "    currently reading is only implemented for"
            print "    data in the Swiss coordinate system"
            print "    x_1 (long name): ", x_1_long_name

        #print "x_1"
        #print x_1.shape
        #print x_1[0] 
        #print x_1[1] - x_1[0] 
        #print x_1[-1]
        #print "y_1"
        #print y_1.shape
        #print y_1[0] 
        #print y_1[1] - y_1[0] 
        #print y_1[-1]
        #quit()

        if hasattr(satscene, 'forecast_time'):
            satscene.forecast_time = append(satscene.forecast_time, ftime)
        else:
            # add first entry 
            satscene.forecast_time = array(ftime)

        # Load specified data sets from netCDF file
        for chn_name in satscene.channels_to_load:

            # Read variable corresponding to channel name
            data = ds.variables[chn_name][:,:,:,:] # use slicing to convert netCDF4 class to masked array         

            # reverse order of the y-axis
            data = data[:, :, ::-1, :]             # dimensions: time, z, y, x

            long_name = ds.variables[chn_name].getncattr("long_name")
            #print "... long_name ", long_name

            FillValue = ds.variables[chn_name].getncattr("_FillValue")
            #print "... FillValue ", FillValue
            data = ma.masked_equal( data, FillValue )

            #print '... '+long_name+".min()/max(): ",data.min(),data.max()
            #print '... '+long_name+".shape: ",   data.shape

            units = ds.variables[chn_name].getncattr("units")
            if units=="1":
                units=None
            #print "... units ", units

            coordinates = ds.variables[chn_name].getncattr("coordinates")

            #from trollimage.image import Image as trollimage
            #img = trollimage(data[0,5,:,:], mode="L", fill_value=[1,1,1]) # [0,0,0] [1,1,1]
            #from trollimage.colormap import rainbow
            #rainbow.set_range(data.min(), data.max())
            #img.colorize(rainbow)
            #img.show()
            #quit()

            if satscene[chn_name].data is None:
                # add first entry 
                satscene[chn_name].data = ma.asarray(data)
                satscene[chn_name].forecast_time = array(ftime)
            else:
                # concatinate two arrays (first entry is time
                satscene[chn_name].data = stack([data[0,:,:,:],satscene[chn_name].data[0,:,:,:]], axis=0)
                satscene[chn_name].forecast_time = append(satscene[chn_name].forecast_time,ftime)
                #print satscene[chn_name].data.shape
                #print satscene[chn_name].forecast_time

            if units != None:
                satscene[chn_name].info['units'] = units

            #satscene.channels.append(chn_name)
