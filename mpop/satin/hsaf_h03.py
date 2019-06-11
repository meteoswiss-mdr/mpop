"""
Reader for EUMETSATs Hydrology SAF (HSAF) h03 product
HSAF website http://hsaf.meteoam.it
h03 product is precipitation rate at the ground 
by GEO(MSG)/Infrared supported by LEO/Microwave
http://hsaf.meteoam.it/precipitation.php?tab=3

After registration the data is available from 
ftp://ftphsaf.meteoam.it/h03

possible accepted formats for this reader are:
* grib as provided by HSAF
* netCDF (grib file converted with cdo)

- Initial version: 
  2015-07-23 Ulrich Hamann (MeteoSwiss)
"""

from mpop import CONFIG_PATH
from ConfigParser import ConfigParser
import os
import numpy.ma as ma
from glob import glob
import datetime

from mpop.utils import get_logger
LOG = get_logger('satin/hsaf_h03')

def find_hsaf_files(time_slot, fullname):

    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, fullname + ".cfg"))
    values = {
        "satname": 'MSG',
        "number": '10',
        "instrument": 'SEVIRI',
        "satellite": fullname
    }

    # HSAF could not decide to take end time 12, 27, 42, 57 or start time 00, 15, 30, 45 or 
    # so we search for both versions here:
    for i_format in [0,1]:
        if i_format==0:
            # time stamp == end of scan time 12min after start, 12, 27, 42, 57, for most of the cases true
            end_time = time_slot + datetime.timedelta(minutes=12)
            filepath    = end_time.strftime(conf.get("seviri-level7", "dir",raw=True))
            filepattern = end_time.strftime(conf.get("seviri-level7", "filename",raw=True)) % values
        else:
            # time stamp == beginning of scan time, for a period between 20170518_0745 and  20170605_1159
            #### datetime_filename_change1 = datetime.datetime(2017, 5, 18,  7, 43, 0)
            #### datetime_filename_change2 = datetime.datetime(2017, 6,  5, 11, 59, 0)
            filepath    = time_slot.strftime(conf.get("seviri-level7", "dir",raw=True))
            filepattern = time_slot.strftime(conf.get("seviri-level7", "filename",raw=True)) % values

        filename_wildcards = os.path.join( filepath, filepattern)
        print "... search for file: ", filename_wildcards
        filenames = glob(str(filename_wildcards))
        
        if len(filenames) == 0:
            print "*** Warning, no HSAF input file found, try to find file with scan start time instead of end time stamp"
        else:
            if len(filenames) > 1:
                print "*** Warning, more than 1 HSAF input data file found"
                for filename in filenames:
                    print "    ", filename
            # found just one (or more) input file, break loop
            break
    
    return filenames
        
        
def load(satscene, **kargs):
    """Reader for EUMETSATs Hydrology SAF (HSAF) h03 product
    h03 product is precipitation rate at the ground 
    by GEO(MSG)/Infrared supported by LEO/Microwave
    http://hsaf.meteoam.it/precipitation.php?tab=3
    """

    filenames = find_hsaf_files(satscene.time_slot, satscene.fullname)
    if len(filenames) == 0:
        print "*** Warning, no HSAF input file found, cannot load."
        #raise IOError("Error, no HSAF input file found, cannot load.")
        LOG.warning("Warning, no HSAF input file found, cannot load.")
        return -1
    else:
        if len(filenames) > 1:
            print "*** Warning, more than 1 HSAF input data file found, choose first one: ", filenames[0]
        filename = filenames[0]
        
    # possible formats: h03_20150513_1557_rom.grb.gz, h03_20150513_1612_rom.grb, h03_20150513_1612_rom.nc
    fileformats = [filename.split(".")[-1] for filename in filenames]
    # try to find grb file 

    if 'grb' in fileformats:
        # read grib 
        data, fill_value, units, long_name = read_h03_grib(filenames[fileformats.index('grb')])
    elif 'nc' in fileformats:
        # read netCDF
        data, fill_value, units, long_name = read_h03_netCDF(filenames[fileformats.index('nc')])
    elif 'gz' in fileformats:
        # unzip 
        from subprocess import call
        from ntpath import basename
        infile = filenames[fileformats.index('gz')]
        outfile = '/tmp/'+basename(infile[:-3])
        print "    unzip ", infile 
        # gunzip -c h03_20150513_1557_rom.grb.gz > h03_20150513_1557_rom.grb
        # call("/bin/gunzip "+ infile                +" 2>&1", shell=True) # dont keep gz file 
        call("/bin/gunzip -c "+ infile+" > "+ outfile  +" 2>&1", shell=True) # keep gz file 
        # check format of gunziped file
        if outfile.split(".")[-1] == 'grb':
            data, fill_value, units, long_name = read_h03_grib(outfile)
        elif outfile.split(".")[-1] == 'nc':
            data, fill_value, units, long_name = read_h03_netCDF(outfile)
        call("rm "+ outfile  +" 2>&1", shell=True) # delete tmporary file 

    if units == "kg m**-2 s**-1" or units == "kg m-2s-1" or units == "kg m-2 s-1":
        data *= 3600 
        units = "kg m-2 h-1"
        
    satscene['h03'] = data
    satscene['h03'].fill_value = fill_value
    satscene['h03'].units      = units
    satscene['h03'].long_name  = long_name
    satscene['h03'].product_name = 'h03'

    # personal communication with help desk
    # Each H03 grib file contains precipitation data of a 900x1900 pixel sub-area of the SEVIRI full disk area (3712x3712 pixels). 
    # The first pixel  of H03 (pixel (1,1)) grib file corresponds to Seviri pixel (1095,85) if the Seviri pixel (1,1) is in the Nort-East. 
    # I can confirm that only the prime satellite is used (position subsatellite longitude 0 degree East).
    # For the future we are thinking to disseminate the h03 outputs already corrected in parallax.

    # conversion of above information to correct AreaDefinition 
    # full_disk = get_area_def("SeviriDiskFull")
    # from mpop.projector import get_area_def
    # import numpy as np
    # np.array(area_def.get_proj_coords(data_slice=(85+900,1095     ))) - 3000.40316582 / 2.
    #    array([-2284807.01076965,  2611850.9558437 ])
    # np.array(area_def.get_proj_coords(data_slice=(85    ,1095+1900))) + 3000.40316582 / 2.
    #                                        array([ 3418959.40744847,  5315214.20824482])
    # or 
    # aex = full_disk.get_area_extent_for_subsets(985,1095,85,2995)

    proj = {'proj': 'geos', 'a': '6378169.0', 'b': '6356583.8', 'h': '35785831.0', 'lon_0': '0.0'}
    aex =      (-2284807.01076965, 2611850.9558437,  3418959.40744847,  5315214.20824482)

    from pyresample.geometry import AreaDefinition
    hsaf_area = AreaDefinition("hsaf", "hsaf", "geos0", proj, 1900, 900, aex)
    
    satscene.area = hsaf_area
    satscene['h03'].area = hsaf_area
    
def read_h03_grib(filename):

    try:
        import pygrib
    except ImportError:
        print "... module pygrib needs to be installed"
        quit()

    # see http://pygrib.googlecode.com/svn/trunk/docs/pygrib-module.html
    print("... read data from %s" % str(filename))

    grbs = pygrib.open(filename)

    #print(grbs)

    #print 'inventory'
    #for grb in grbs:
    #    print(grb)
    #print 'end inventory'

    long_name  = 'Instantaneous rain rate'
    units      = 'kg m**-2 s**-1'
    _FillValue = 0.0 
    grb = grbs.select(name=long_name)[0]
    # print(grb)
    data = ma.asarray(grb.values)
    
    # mask data for no rain and "no cloud" marker
    data = ma.masked_where(data ==    0.0, data)
    data = ma.masked_where(data == 9999.0, data)

    print '    fill_value: ', 0 
    print '    units:      ', units
    print '    long_name:  ', long_name 
    print '    datatype:   ', type(data)
    print '    shape:      ', data.shape
    print '    min/max:    ', data.min(), data.max()

    return data, _FillValue, units, long_name 


def read_h03_netCDF(filename):

    try:
        from netCDF4 import Dataset
    except ImportError:
        print "... module netCDF4 needs to be installed"
        quit()

    print("... read data from %s" % str(filename))

    # Load data from netCDF file
    # see also http://netcdf4-python.googlecode.com/svn/trunk/docs/netCDF4-module.html
    ds = Dataset(filename, 'r')

    if 'irrate' in ds.variables:
        print '    found variable irrate'
        var_name='irrate'            # converted with: cdo -f nc4c -z zip copy infile outfile
    elif 'IRRATE_P30_GSV0' in ds.variables:
        print '    found variable IRRATE_P30_GSV0'
        var_name='IRRATE_P30_GSV0'   # converted with: ncl_convert2nc h03_20150529_0827_rom.grb -e grb -nc4c -cl 9
        #variables:
        #    float IRRATE_P30_GSV0(ygrid_0, xgrid_0) ;
        #            IRRATE_P30_GSV0:initial_time = "05/29/2015 (08:27)" ;
        #            IRRATE_P30_GSV0:parameter_template_discipline_category_number = 30, 3, 1, 1 ;
        #            IRRATE_P30_GSV0:parameter_discipline_and_category = "Space products, Quantitative products" ;
        #            IRRATE_P30_GSV0:grid_type = "Space view perspective or orthographic" ;
        #            IRRATE_P30_GSV0:_FillValue = 1.e+20f ;
        #            IRRATE_P30_GSV0:units = "kg m-2s-1" ;
        #            IRRATE_P30_GSV0:long_name = "Instantaneous rain rate" ;
        #            IRRATE_P30_GSV0:production_status = "Operational test products" ;
        #            IRRATE_P30_GSV0:center = "Rome (RSMC)" ;
        print '*** Error, does not work for unknown reason'
        print '    data.mask = (data == _FillValue) | (data == 0.0) produce error'
        quit()

    #print type(ds.variables[var_name])
    #print dir(ds.variables[var_name])
  
    _FillValue = ds.variables[var_name]._FillValue
    # or fill_value = ds.variables[var_name].getncattr('_FillValue')
    units      = ds.variables[var_name].units
    long_name  = ds.variables[var_name].long_name 

    # Read variable corresponding to channel name
    data = ma.asarray(ds.variables[var_name])

    print '    fill_value: ', ds.variables[var_name]._FillValue 
    print '    units:      ', ds.variables[var_name].units
    print '    long_name:  ', ds.variables[var_name].long_name 
    print '    datatype:   ', ds.variables[var_name].datatype
    print '    shape:      ', data.shape
    print '    min/max:    ', data.min(), data.max()

    if len(data.shape) == 3:
        if data.shape[0] == 1:
            print "   reduce to 2 dimensions (skip time dimension)"
            data = ma.asarray(ds.variables[var_name][0,:,:])
        else:
            print "*** Error, unknown netCDF file format in h03_nc.py"
            print "    probably more time steps in one file (not implemented yet)"
            quit()

    # mask data for no rain, fill_value and no cloud marker
    data = ma.masked_where(data ==        0.0, data)
    data = ma.masked_where(data == _FillValue, data)
    data = ma.masked_where(data ==     9999.0, data)

    return data, _FillValue, units, long_name 
