import Image
import logging
import glob
import os
from ConfigParser import ConfigParser
import glob
import numpy as np
import numpy.ma as ma
from mpop import CONFIG_PATH

import pyresample
import datetime
from time import gmtime

from os.path import dirname, exists
from os import makedirs

import logging
#LOG = logging(__name__)



def load(satscene, *args, **kwargs):
    """Loads the *channels* into the satellite *scene*.
    """

    # Dataset information
    #
    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))

    values = {"orbit": satscene.orbit,
              "satname": satscene.satname,
              "number": satscene.number,
              "instrument": satscene.instrument_name,
              }

    # print values
    if "area" in kwargs:
        projectionName = kwargs["area"]
    else:
        projectionName = conf.get("Overshooting_Tops-level2", "projection")
        #projectionName = "ccs4"
    print "... project lightning data to ", projectionName
    #scale = conf.get("radar-1", "scale")
    #print '... read scale from: ', scale

    timeslot_end_obs = satscene.time_slot + datetime.timedelta(minutes=4)

    filename = os.path.join(
        timeslot_end_obs.strftime(conf.get("Overshooting_Tops-level2", "dir")),
        timeslot_end_obs.strftime(conf.get("Overshooting_Tops-level2", "filename", raw=True)) % values)
   
    # Load data from txt file 
    print "... read data from", str(filename)
    if len(glob.glob(str(filename))) == 0:
        "*** WARNING, no file "+str(filename)+" found!"
        filename=""
    elif len(glob.glob(str(filename))) > 1:
        "*** WARNING, more than one file "+str(filename)+" found!"
        filename = glob.glob(str(filename))[0]
    else: 
        filename = glob.glob(str(filename))[0]

    from netCDF4 import Dataset
    fh = Dataset(filename, mode='r')
        
    print "... channels to load ", satscene.channels_to_load
    for chn_name in satscene.channels_to_load:
        # Read variable corresponding to channel name
        print "... channel to read ", chn_name
        if chn_name == 'longitude':
            satscene[chn_name] = fh.variables['longitude'][:]
        if chn_name == 'latitude':
            satscene[chn_name] = fh.variables['latitude'][:]
        elif chn_name == 'ir_brightness_temperature':
            satscene[chn_name] = fh.variables['ir_brightness_temperature'][0,:,:]
        elif chn_name == 'ot_potential_temperature':
            satscene[chn_name] = fh.variables['ot_potential_temperature'][0,:,:]
        elif chn_name == 'ot_height':
            satscene[chn_name] = fh.variables['ot_height'][0,:,:]
        elif chn_name == 'ot_pressure':
            satscene[chn_name] = fh.variables['ot_pressure'][0,:,:]
        elif chn_name == 'ot_anvilmean_brightness_temperature_difference':
            satscene[chn_name] = fh.variables['ot_anvilmean_brightness_temperature_difference'][0,:,:]
        elif chn_name == 'position_parallax_corrected':
            lons_pc = fh.variables['parallax_correction_longitude'][0,:,:]
            lats_pc = fh.variables['parallax_correction_latitude'][0,:,:]
            from mpop.projector import get_area_def
            obj_area = get_area_def(projectionName)
            #dens = lon_lat_to_dens(lons, lats, obj_area)
            dens_pc = lon_lat_to_dens(lons_pc.flatten(), lats_pc.flatten(), obj_area)
            satscene[chn_name] = dens_pc
            # get projection 
            projection = pyresample.utils.load_area(os.path.join(CONFIG_PATH, "areas.def"), projectionName)
            #print projection
            satscene.area = projection
            

# ---------------------------------------------------------------------------------

def lon_lat_to_dens(lons, lats, obj_area):

    ni = obj_area.y_size # 640 for ccs4
    nj = obj_area.x_size # 710 for ccs4
   
    dens = ma.asarray(np.zeros(shape=(ni,nj)))

    (jCH, iCH) = obj_area.get_xy_from_lonlat(lons, lats, outside_error=False)       # for Bern (345, 280)
    # returns an int, if only one element lon and lat data is passed
    if isinstance(jCH, int):
        jCH = [jCH]
    if isinstance(iCH, int):
        iCH = [iCH]

    for i, j in zip(iCH, jCH):
        #print "   ", str(t_light), " ", str(time_slot), " ", str(time_slot-t_light), "   ", ltype, "          ", curr
        if 0<i and i<ni and 0<j and j<nj:
            dens[i,j]+=1

    return dens

# ---------------------------------------------------------------------------------


def readLightning(file, NEAR_REAL_TIME, time_slot, dt=5, area='ccs4'):
    """ input 
        file           [string]   input file name
        NEAR_REAL_TIME [logical]  True  -> data from near real time folder
                                  False -> data from archive
        time_slot      [datetime] desired time of lightning analysis 
        dt             [float]    period in min a lightning is counted for
        form           [string]   'square' -> lightning is counted in a square
                                  'circle' -> lightning is counted in a circle
    """

    from numpy import genfromtxt, floor, ceil, zeros, nonzero, rot90, flipud
    from ApproxSwissProj import ApproxSwissProj 
    import numpy.ma as ma
    from os import stat

    from mpop.projector import get_area_def
    obj_area = get_area_def(area)

    ni = obj_area.y_size # 640 for ccs4
    nj = obj_area.x_size # 710 for ccs4

    # empty masked array
    dens     = ma.asarray(zeros(shape=(ni,nj)))
    densIC   = ma.asarray(zeros(shape=(ni,nj)))
    densCG   = ma.asarray(zeros(shape=(ni,nj)))
    curr_abs = ma.asarray(zeros(shape=(ni,nj)))
    curr_neg = ma.asarray(zeros(shape=(ni,nj)))
    curr_pos = ma.asarray(zeros(shape=(ni,nj)))

    # read data from file
    if file != "" and stat(file).st_size != 0:
        
        if NEAR_REAL_TIME:
            # lat        lon    time           intra  current
            # 43.2540    9.9290 20150201033521 IC     35.2
            #                   12345678901234
            data = genfromtxt(file, dtype=None, names=('lat', 'lon', 'time', 'intra', 'current'))  # automatic specification of the format
            if data['time'].size == 1:   # here I get an error for the implecit loop: ??? TypeError: iteration over a 0-d array ???
                years  = [str(data['time'])[ 0: 4]]
                months = [str(data['time'])[ 4: 6]]
                days   = [str(data['time'])[ 6: 8]]
                hours  = [str(data['time'])[ 8:10]]
                mins   = [str(data['time'])[10:12]]
                secs   = [str(data['time'])[12:14]]
                intra   = [data['intra']]   # just a string, instead of an one element string array
                current = [data['current']] # just a string, instead of an one element string array

            else:
                years  = [str(date)[ 0: 4] for date in data['time']]
                months = [str(date)[ 4: 6] for date in data['time']]
                days   = [str(date)[ 6: 8] for date in data['time']]
                hours  = [str(date)[ 8:10] for date in data['time']]
                mins   = [str(date)[10:12] for date in data['time']]
                secs   = [str(date)[12:14] for date in data['time']]
                intra   = data['intra']   
                current = data['current']
            #intra_cloud = zeros(len(data['type']), dtype=bool) # first initialize all with False == cloud to ground lightning
            #intra_cloud[nonzero(data['type'])] = True  
        else: 
            # date                     |  lon   |  lat   |  curr|nS|mode|intra|    Ax |   Ki2|   Exc|      Incl|    Arc|  empty|  empty|  empty|  empty|
            # 01.09.2013 12:19:39.0 UTC|  12.824|  46.506| -16.0| 1|   6|    0|    0.1|   0.4|   1.0|       145|    1.0|    0.0|    0.0|    0.0|    0.0|
            data = genfromtxt(file, dtype=None, names=('date', 'lon', 'lat','current','nS','mode','intra','Ax','Ki2','Exc','Incl','Arc','d1','d2','d3','d4','d5'), delimiter="|") # automatic specification of the format
            days   = [str(date)[ 0: 2] for date in data['date']]
            months = [str(date)[ 3: 5] for date in data['date']]
            years  = [str(date)[ 6:10] for date in data['date']]
            hours  = [str(date)[11:13] for date in data['date']]
            mins   = [str(date)[14:16] for date in data['date']]
            secs   = [str(date)[17:21] for date in data['date']]
            intra   = data['intra']   
            current = data['current']

        #if area=='ccs4':
        if False:
            # define class for reprojection to ccs4 grid 
            ccs4 = ApproxSwissProj()

            # for Bern 
            # data['lat'] = 46.9479
            # data['lon'] = 7.4446
 
            # position in km (x->[-160:480], y->[255,965]
            xCH = ccs4.WGStoCHx(data['lat'],data['lon']) # and x in the vertical         # for Bern 199646.19
            yCH = ccs4.WGStoCHy(data['lat'],data['lon']) # y is in horizontal direction  # for Bern 600454.29

            ## convert to index (i->[0,640], j->[0,710]
            #iCH = floor(yCH/1000.)-255  # y is in horizontal direction
            #jCH = 710-(480-floor(xCH/1000.))  # and x in the vertical ## 
            # convert to index (i->[0,640], j->[0,710]
            iCH = (480-ceil(xCH/1000.))  ## 710- ## x/i is in vertical direction         # for Bern 280.0
            jCH = floor(yCH/1000.)-255    # y/j is in horizontal dir.                    # for Bern 345.0
        else:
            print "... new experimental reprojection"

            if True:
                (jCH, iCH) = obj_area.get_xy_from_lonlat(data['lon'], data['lat'], outside_error=False)       # for Bern (345, 280)
                # returns an int, if only one element lon and lat data is passed
                if isinstance(jCH, int):
                    jCH = [jCH]
                if isinstance(iCH, int):
                    iCH = [iCH]
            else:
                from pyproj import Proj
                p = Proj(obj_area.proj4_string)
                (yCH, xCH) = p(data['lon'], data['lat'])                                 # for Bern (600381.87, 199499.18)
                y_offset = obj_area.area_extent[3]/obj_area.pixel_size_y                 #                 for ccs4 480.0
                iCH = y_offset - ceil(xCH/1000.)                                         # for Bern 280.0
                x_offset = obj_area.area_extent[0]/obj_area.pixel_size_x                 #                 for ccs4 255.0 
                jCH = floor(yCH/1000.) - x_offset                                        # for Bern 345.0

        # considered time period 
        dtime = datetime.timedelta(minutes=dt)

        print "--- lightning time        desired time slot     dtime    intra cloud     current[kA]"
        for i, j, YYYY, MM, DD, hh, mm, ss, ltype, curr in zip(iCH, jCH, years, months, days, hours, mins, secs, intra, current):
            t_light = datetime.datetime(int(YYYY), int(MM), int(DD), int(hh), int(mm), int(float(ss)) ) 
            if ( t_light < time_slot and time_slot-t_light < dtime ):
                #print "   ", str(t_light), " ", str(time_slot), " ", str(time_slot-t_light), "   ", ltype, "          ", curr
                if 0<i and i<ni and 0<j and j<nj:
                    dens[i,j]+=1
                    curr_abs[i,j]+=abs(curr)
                    if curr < 0:
                        curr_neg[i,j]+=curr
                    else:
                        curr_pos[i,j]+=curr
                    if ltype == 'IC' or ltype == 1:
                        #print "... only intra clouds lightnings", YYYY, MM, DD, hh, mm, ltype
                        densIC[i,j]+=1  
                    if ltype=='CG' or ltype==0:
                        #print "... only clouds2ground lightnings", YYYY, MM, DD, hh, mm, ltype
                        densCG[i,j]+=1
            if time_slot + dtime < t_light :
                break
    else:
        print '*** Warning, empty lightning input file '

    print '... dens ndim (lightning.py): ', dens.ndim
    print '... shape (lightning.py): ', dens.shape
    print '... min/max (lightning.py): ', dens.min(), dens.max()

    return dens, densIC, densCG, curr_abs, curr_neg, curr_pos

# ---------------------------------------------------------------------------------

def add_lightning(prop, i, j, dx, form):
    """ add one lightning as data field to the ccs4 array
        add_lightning(dens, i, j, dx, form)
        form           [string]   'square' -> lightning is counted in a square
                                  'circle' -> lightning is counted in a circle
    """

    # print "add lightning at", i, j
    if form == 'square':
        # mark a square 
        prop[i-dx:i+dx,j-dx:j+dx]+=1
    elif form == 'circle':
        from numpy import arange, meshgrid
        # mark a circle 
        rr=dx  # radius
        x = arange(0, 710)  # 709-
        y = arange(0, 640)
        X, Y = meshgrid(x, y)
        #print X
        #print Y        
        #print X.shape
        #print Y.shape
        #print lightnings.shape
        interior = ((X-j)**2 + (Y-i)**2) < rr**2
        prop+=1*interior
    else:
        print '*** ERROR in add_lightning (lightning.py)'
        print '    unknown form ', form

# ---------------------------------------------------------------------------------

def unfold_lightning(prop, dx, form):
    """ unfold one lightning in such a way that it effects dx km
        form           [string]   'square' -> lightning is counted in a square
                                  'circle' -> lightning is counted in a circle
    """

    from scipy.ndimage.filters import convolve
    from numpy import zeros, ones
    import numpy.ma as ma


    # print "unfold_lightning at", i, j
    if form == 'square':
        if dx != 1:
            print "... counting lightning per square with edge lengths of ", dx, " km"
            # mark a square 
            k = ones([2*dx-1,2*dx-1])
            return ma.asarray( convolve(prop, k, mode="constant", cval=0) ) 

    elif form == 'circle':
        if dx != 1:
            print "... counting lightning inside a circle with a radius of ", dx, " km"
            from numpy import arange, meshgrid
            # mark a circle 
            x = arange(0, 2*dx-1) 
            X, Y = meshgrid(x, x)
            k = 1 * ((X-(dx-1))**2 + (Y-(dx-1))**2) < dx**2
            return ma.asarray( convolve(prop, k, mode="constant", cval=0) ) 

    elif form == 'gauss':
        from scipy import ndimage
        print "unfold_lightning (lighting.py) min/max = ", prop.min(), prop.max()
        print '... apply gauss filter with half width of (0.5*', dx, ") km"
            # tuned that it looks similar to the circle with the same dx
        data = ma.asarray( ndimage.gaussian_filter(prop, 0.5*dx) )
        print "unfold_lightning (lighting.py) min/max = ", data.min(), data.max()
        return ma.masked_less(data, 0.001)

    elif form == 'gaussgauss':
        from scipy import ndimage
        print "unfold_lightning (lighting.py) min/max = ", prop.min(), prop.max()
        print '... apply gauss filter with half width of (0.5*', dx, ") km"
            # tuned that it looks similar to the circle with the same dx
        prop = ndimage.gaussian_filter(prop, 0.5*dx)
        data = ma.asarray( ndimage.gaussian_filter(prop, 0.5*dx) )
        print "unfold_lightning (lighting.py) min/max = ", data.min(), data.max()
        return ma.masked_less(data, 0.001)


    else:
        print '*** ERROR in unfold_lightning (lightning.py)'
        print '    unknown form ', form

