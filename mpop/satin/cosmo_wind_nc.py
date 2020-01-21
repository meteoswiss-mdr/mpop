"""Loader for MSG, netcdf format.
"""
from __future__ import division
from __future__ import print_function

from ConfigParser import ConfigParser
from mpop import CONFIG_PATH
import os
from netCDF4 import Dataset
import numpy.ma as ma
from glob import glob
from mpop.projector import get_area_def
from datetime import timedelta, datetime
from numpy import logical_or, flipud, array, stack, append, where, array

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
        print("... %s: %s" % (key, kwargs[key]))
    else:
        forecast_time = [0]

    var_to_load=[]
    if "pressure_levels" in kwargs:
        key = "pressure_levels"
        pressure_levels = array(kwargs[key])
        print("... %s: %s" % (key, kwargs[key]))
    else:
        pressure_levels=[]
        
    for chn_name in satscene.channels_to_load:
        if chn_name[0]=="U":            
            if "U" not in var_to_load:
                var_to_load.append("U")
        if chn_name[0]=="V":            
            if "V" not in var_to_load:
                var_to_load.append("V")
        if len(chn_name)==8:
            if int(chn_name[2:5]) not in pressure_levels:
                pressure_levels.append(int(chn_name[2:5]))

    # convert from hPa to Pa
    if len(pressure_levels)==0:
        print("ERROR, no pressure levels specified when reading COSMO wind")
        quit()
    else:
        print("read ", pressure_levels)

    pressure_levels = array(pressure_levels)
    pressure_levels *= 100
    
    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))
    values = {"satname": satscene.satname,
    "number": satscene.number,
    "instrument": satscene.instrument_name,
    }

    # get COSMO model start time
    hour_run = satscene.time_slot.hour //3 * 3 
    t_run = datetime(satscene.time_slot.year, satscene.time_slot.month, satscene.time_slot.day, hour_run, 0)
    print("COSMO model start", str(t_run))
    
    dt = satscene.time_slot - t_run
    hour_forecast1 = "%02d" % int (dt.total_seconds() / 3600) # using integer devision 
    hour_forecast2 = "%02d" % int (dt.total_seconds() / 3600 +1)  # using integer devision 

    print("hour_forecast1, hour_forecast2", hour_forecast1, hour_forecast2)

    t0 = satscene.time_slot
    print("time to show:", t0)
    ftime = int(satscene.time_slot.strftime("%H")) % 3
    model_starthour = int(satscene.time_slot.strftime("%H")) - ftime
    model_starttime = datetime(t0.year, t0.month, t0.day, model_starthour, 0, 0) 
    print("available calc time (until now): ", datetime.now()-model_starttime)
    # check if COSMO model already produced output, if not take run from 3 hours before
    if datetime.now()-model_starttime < timedelta(minutes=80):
        print("... COSMO results from model start ",model_starttime ," not yet ready,. take run three hours before")
        model_starttime -= timedelta(minutes=180)
        ftime += 3
    print(model_starttime)
    forecast_time = [ftime]
    
    for ftime in forecast_time:

        ftime_str = "%02d" % ftime # in hours
        values["forecast_time"] = ftime_str
        print('... forecast time: '+ftime_str)

        filename = os.path.join( model_starttime.strftime(conf.get("cosmo-level3", "dir", raw=True)),
                                 model_starttime.strftime(conf.get("cosmo-level3", "filename", raw=True)) % values )

        print("... search for file: ", filename)
        filenames=glob(str(filename))
        if len(filenames) == 0:
            print("*** Error, no file found")
            quit()
        elif len(filenames) == 1:
            filename = filenames[0]
            # if bz2 file then bunzip2 the file
            if filename.split(".")[-1]=="bz2":
                from subprocess import call
                from ntpath import basename
                #infile = filenames[fileformats.index('bz2')]
                call("cp "+ filename+" /tmp 2>&1", shell=True)
                tmpfile = '/tmp/'+basename(filename)
                print("... bunzip2 "+tmpfile)
                call("/bin/bunzip2 "+ tmpfile+" 2>&1", shell=True)
                print("... remove "+tmpfile)
                call("rm "+ tmpfile+" 2>&1", shell=True)
                filename = tmpfile[:-4]
                copy_file=True
        elif len(filenames) > 1:
            print("*** Warning, more than 1 datafile found: ", filenames)
            for this_file in filenames:
                if this_file.split(".")[-1] == 'nc':
                    filename=this_file
                    print("*** Choose first nc file: ", filename)
                    break
        
        print(("... read data from %s" % str(filename)))
        
        # Load data from netCDF file
        ds = Dataset(filename, 'r')

        # read the dimensions
        x_1 = ds.variables['x_1']
        y_1 = ds.variables['y_1']
        z_1 = ds.variables['z_1']
        satscene.x_1 = x_1[:]
        satscene.y_1 = y_1[:]
        satscene.z_1 = z_1[:]
        satscene.x_1_long_name = ds.variables['x_1'].getncattr("long_name")
        satscene.y_1_long_name = ds.variables['y_1'].getncattr("long_name")
        satscene.z_1_long_name = ds.variables['z_1'].getncattr("long_name")

        #lon_1 = ds.variables['lon_1']
        #lat_1 = ds.variables['lat_1']

        #print x_1[:], len(x_1)
        #print y_1[:], len(y_1)
        print(z_1[:])

        if hasattr(satscene, 'forecast_time'):
            satscene.forecast_time = append(satscene.forecast_time, ftime)
        else:
            # add first entry 
            satscene.forecast_time = array(ftime)
            
        # Load specified data sets from netCDF file
        for chn_name in var_to_load:

            print("... read:", chn_name)
            
            # Read variable corresponding to channel name
            data = ds.variables[chn_name][:,:,:,:] # use slicing to convert netCDF4 class to masked array         

            # reverse order of the y-axis
            data = data[:, :, ::-1, :]             # dimensions: time, z, y, x

            long_name = ds.variables[chn_name].getncattr("long_name")
            #print "... long_name ", long_name

            FillValue = ds.variables[chn_name].getncattr("_FillValue")
            #print "... FillValue ", FillValue
            data = ma.masked_equal( data, FillValue )

            print('... '+long_name+".min()/max(): ",data.min(),data.max())
            print('... '+long_name+".shape: ",   data.shape)

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

            #satscene.channels.append(chn_name)

            for p_level in pressure_levels:

                print('... read wind for pressure level:', p_level)

                if p_level in z_1:
                    z_index = where(z_1[:]==p_level)
                else:
                    print("*** ERROR while reading COSMO wind")
                    print("    unknown pressure level", p_level)
                    quit()

                p_str= '%d' % (p_level/100)
                p_str = p_str.strip() +'hPa'
                print(p_str)
                var_name=chn_name+'-'+p_str
                #var_name=chn_name
                
                satscene[var_name].data = ma.asarray(data[0,int(z_index[0]),:,:])
                satscene[var_name].forecast_time = array(ftime)

                #if chn_name not in satscene:
                #    print "add first data entry"
                #    print data.shape, data[0,0,:,:].shape
                #    satscene[chn_name].data = ma.asarray(data[0,0,:,:])
                #    satscene[chn_name].forecast_time = array(ftime)
                #else:
                #    print "concatinate two arrays (first dimension is time)"
                #    satscene[chn_name].data = stack([data[0,:,:,:],satscene[chn_name].data[0,:,:,:]], axis=0)
                #    satscene[chn_name].forecast_time = append(satscene[chn_name].forecast_time,ftime)
                #    #print satscene[chn_name].data.shape
                #    #print satscene[chn_name].forecast_time

                if units != None:
                    satscene[var_name].info['units'] = units

                if satscene.x_1_long_name != "x-coordinate in Swiss coordinate system" :
                    print("*** Error in mpop/satin/cosmo2_wind_nc.py")
                    print("    currently reading is only implemented for")
                    print("    data in the Swiss coordinate system")
                    print("    x_1 (long name): ", x_1_long_name)
                    quit()

                # !!! BAD: THIS INFORMATION SHOULD BE SAVED IN THE FILE !!!
                if len(x_1)==1000 and len(y_1)==700:
                    satscene[var_name].area = get_area_def("swissXXL")
                elif len(x_1)==710 and len(y_1)==640:
                    satscene[var_name].area = get_area_def("ccs4")
                else:
                    print("*** ERROR")
                    print("unknown area of COSMO input")
                    quit()
