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

long_names = {0.6:'VIS006', 0.8:'VIS008', 1.6:'IR_016', 3.9:'IR_039', 6.2:'WV_062', 7.3:'WV_073',\
              8.7:'IR_087', 9.7:'IR_097', 10.8:'IR_108', 12.0:'IR_120', 13.4:'IR_134',\
              'VIS006':'VIS006', 'VIS008':'VIS008', 'IR_016':'IR_016', 'IR_039':'IR_039', 'WV_062':'WV_062', 'WV_073':'WV_073',\
              'IR_087':'IR_087', 'IR_097':'IR_097', 'IR_108':'IR_108', 'IR_120':'IR_120', 'IR_134':'IR_134', 'HRV':'HRV'}
units      = {0.6:'percent', 0.8:'percent', 1.6:'percent', 3.9:'K', 6.2:'K', 7.3:'K',\
              8.7:'K', 9.7:'K', 10.8:'K', 12.0:'K', 13.4:'K',\
              'VIS006':'percent', 'VIS008':'percent', 'IR_016':'percent', 'IR_039':'K', 'WV_062':'K', 'WV_073':'K',\
              'IR_087':'K', 'IR_097':'K', 'IR_108':'K', 'IR_120':'K', 'IR_134':'K', 'HRV':'percent'}

def load(satscene, **kwargs):
    """Load MSG SEVIRI radiances from netCDF file.
    """
    
    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))

    if "reader_level" in kwargs.keys():
        if kwargs["reader_level"] is not None and kwargs["reader_level"] != "":
            reader_level = kwargs["reader_level"]
        else:
            reader_level = "seviri-level8"
    else:
        reader_level = "seviri-level8"
    print "... use reader_level:", reader_level
                
    # assume this is reader_level seviri-level8
    filename = os.path.join( satscene.time_slot.strftime(conf.get(reader_level, "dir", raw=True)),
                             satscene.time_slot.strftime(conf.get(reader_level, "filename", raw=True)))

    print "... search for file: ", filename
    filenames=glob(str(filename))
    if len(filenames) == 0:
        print "*** Error, no file found"
        return -1
    elif len(filenames) > 1:
        print "*** Warning, more than 1 datafile found: ", filenames 
    filename = filenames[0]
    print("... read data from %s" % str(filename))

    # Load data from netCDF file
    ncfile = Dataset(filename,'r')
        
    # !!! BAD: THIS INFORMATION SHOULD BE READ FROM THE FILE !!!
    #       int64 grid_mapping_0 ;
    #            grid_mapping_0:ellipsoid = "bessel" ;
    #            grid_mapping_0:false_easting = 600000LL ;
    #            grid_mapping_0:false_northing = 200000LL ;
    #            grid_mapping_0:grid_mapping_name = "swiss_oblique_mercator" ;
    #            grid_mapping_0:latitude_of_projection_origin = 46.9524055555556 ;
    #            grid_mapping_0:longitude_of_projection_origin = 7.43958333333333 ;
    #            grid_mapping_0:scale_factor_at_projection_origin = 1LL ;
    #            grid_mapping_0:units = "m" ;
    # !!! CURRENTLY: area is set fixed to ccs4     
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
