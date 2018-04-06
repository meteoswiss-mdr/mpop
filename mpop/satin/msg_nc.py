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

long_names = {3.9:'IR_039',9.7:'IR_097', 7.3:'WV_073', 6.2:'WV_062', 10.8:'IR_108','IR_039':'IR_039', 'IR_097':'IR_097', 'WV_073':'WV_073', 'WV_062':'WV_062', 'IR_108':'IR_108', 'IR_120':'IR_120', 'HRV':'HRV'}
units      = {9.7: 'K', 7.3: 'K', 6.2: 'K', 10.8: 'K', 'IR_097': 'K', 'WV_073': 'K', 'WV_062': 'K', 'IR_108': 'K', 'IR_120':'K'}

def load(satscene, **kargs):
    """Load MSG SEVIRI radiances from netCDF file.
    """
    
    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))

    # assume this is reader_level seviri-level8
    filename = os.path.join( satscene.time_slot.strftime(conf.get("seviri-level8", "dir", raw=True)),
                             satscene.time_slot.strftime(conf.get("seviri-level8", "filename", raw=True)))

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
        if chn_name in long_names.keys():
        
            # Read variable corresponding to channel name
            data = ncfile.variables[chn_name][:,:,0] # attention [:,:] or [:] is really necessary
            #print type(data)                         # it converts 'NetCDFVariable' to 'numpy.ndarray'
            #print data.ndim
            #print data.shape

            data_m = ma.masked_array(data) # convert numpy array to masked array 
            data_m.mask = (data == 999.0)  # create mask 

            # with this rotation, we can use same indices as we do in the EUMETSAT order mask 
            # data_m.mask[3221:3371, 1712:2060] = 1

            #import matplotlib.pyplot as plt
            #plt.imshow(data_m)
            #plt.colorbar()
            #plt.show()

            satscene[chn_name] = data_m   #  
            satscene[chn_name].fill_value = 999.0
            #satscene[chn_name].info['units'] = units[chn_name]
            satscene[chn_name].info['units'] = 'K'
            satscene[chn_name].long_name  = long_names[chn_name]
            satscene[chn_name].area = area_def

    # close the file.
    ncfile.close()        
