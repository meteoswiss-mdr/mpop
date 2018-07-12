"""Loader for MSG, netcdf format.
"""
from ConfigParser import ConfigParser
from mpop import CONFIG_PATH
import os

from netCDF4 import Dataset

import numpy.ma as ma
from glob import glob
from mpop.projector import get_area_def
from numpy import flipud, ndarray
from datetime import datetime

cosmo_vars_2d=["lon_1", "lat_1",\
               "TWATER", "tropopause_height", "tropopause_temperature", "tropopause_pressure", \
               "FF_10M", "VMAX_10M", "CAPE_MU", "CAPE_ML", "CIN_MU", "CIN_ML", \
               "SLI", "LCL_ML", "LFC_ML", "T_2M", "TD_2M", "GLOB", "PS", \
               "PMSL", "PMSLr", "HZEROCL", "WSHEAR_0-3km", "WSHEAR_0-6km", "SYNMSG_BT_CL_IR10.8"]
cosmo_vars_3d=["POT_VORTIC", "THETAE", "MCONV", "geopotential_height",\
               "T_SO", "RELHUM"]



def load(satscene, **kwargs):
    """Load COSMO forecasts from netCDF file.
    """
    
    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))

    t0=satscene.time_slot
    print t0
    hour2show=int(satscene.time_slot.strftime("%H"))
    print hour2show, type(hour2show)
    ftime = int(hour2show) % 3
    print ftime
    model_starthour = hour2show - ftime
    model_starttime = datetime(t0.year, t0.month, t0.day, model_starthour, t0.minute, 0) 
    print model_starttime
    
    ftime_str = "%02d" % ftime # in hours
    values={}
    values["forecast_time"] = ftime_str
    
    # assume this is reader_level cosmo-level2
    filename = os.path.join( model_starttime.strftime(conf.get("cosmo-level2", "dir", raw=True)),
                             model_starttime.strftime(conf.get("cosmo-level2", "filename", raw=True)) % values)

    print "... search for file: ", filename
    filenames=glob(str(filename))
    if len(filenames) == 0:
        print "*** Error, no file found"
        quit()
    elif len(filenames) == 1:
        filename = filenames[0]
        #if filename.split(".")[-1]=="bz2":
        # from subprocess import call
        # from ntpath import basename
        # infile = filenames[fileformats.index('bz2')]
        # call("cp "+ infile+" /tmp 2>&1", shell=True)
        # tmpfile = '/tmp/'+basename(infile)
        # print "... bunzip2 "+tmpfile
        # call("/bin/bunzip2 "+ tmpfile+" 2>&1", shell=True)
        # print "... remove "+tmpfile
        # call("rm "+ tmpfile+" 2>&1", shell=True)
        # filename = tmpfile[:-4]
        # copy_file=True
    elif len(filenames) > 1:
        print "*** Warning, more than 1 datafile found: ", filenames
        for this_file in filenames:
            if this_file.split(".")[-1] == 'nc':
                filename=this_file
                print "*** Choose first nc file: ", filename
                break
        
    print("... read data from %s" % str(filename))

    # Load data from netCDF file
    ncfile = Dataset(filename,'r')
        
    # area is set fixed to ccs4 
    # !!! BAD: THIS INFORMATION SHOULD BE SAVED IN THE FILE !!!
    area = 'ccs4'
    area_def = get_area_def(area)
    print "... set area to ", area 
    #satscene.area = area_def
        
    for chn_name in satscene.channels_to_load:
        #print '    load ', chn_name
        
        # only read variables which are in the netCDF file
        if chn_name in cosmo_vars_2d+cosmo_vars_3d:

            print ncfile.variables.keys()
            print type(ncfile.variables[chn_name])
            print ncfile.variables[chn_name].shape
            
            # Read variable corresponding to channel name
            data = ncfile.variables[chn_name][0,:,:] # attention [:,:] or [:] is really necessary
            print type(data)                         # it converts 'NetCDFVariable' to 'numpy.ndarray'
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

            satscene[chn_name] = flipud(data_m)   #  
            #satscene[chn_name].fill_value    = ncfile.variables[chn_name].attributes()["_FillValue"] 
            #satscene[chn_name].info['units'] = ncfile.variables[chn_name].attributes()["units"]
            #satscene[chn_name].long_name     = ncfile.variables[chn_name].attributes()["long_name"]
            satscene[chn_name].area = area_def

    # close the file.
    ncfile.close()        
