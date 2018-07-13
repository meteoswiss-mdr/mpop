from PIL import Image
#import Image
import logging
import glob
import os
from ConfigParser import ConfigParser
import glob
import numpy as np
import numpy.ma as ma
from mpop import CONFIG_PATH

import pyresample

import logging
#LOG = logging(__name__)

verbose=True

def load(satscene, *args, **kwargs):
    """Loads the *channels* into the satellite *scene*.
    """
    #
    # Dataset information
    #
    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))

    values = {"orbit": satscene.orbit,
              "satname": satscene.satname,
              "number": satscene.number,
              "instrument": satscene.instrument_name,
              "satellite": satscene.fullname
              }
    # print values
    projectionName = conf.get("radar-level2", "projection")

    for chn_name in satscene.channels_to_load:

        for i in xrange(9):
            radar_product='radar-{0:1d}'.format(i+1)
            prod_name = conf.get(radar_product, "name")
            print prod_name, chn_name
            if prod_name.replace("'", "").replace('"', '') == chn_name.replace("'", "").replace('"', ''):
                values["product"]=conf.get(radar_product, "product_name").replace("'", "").replace('"', '')
                if verbose:
                    print "... radar product to read:", chn_name, ' (', values["product"], ')'
                scale = conf.get(radar_product, "scale")
                if verbose:
                    print '... read scale from: ', scale
                units = conf.get(radar_product, "units").replace("'", "").replace('"', '')
                break

        filename = os.path.join(
            satscene.time_slot.strftime(conf.get("radar-level2", "dir", raw=True)) % values,
            satscene.time_slot.strftime(conf.get("radar-level2", "filename", raw=True)) % values)

        # Load data from file
        # filename = "/home/lom/users/cll/pytroll/radar/RZC143081530F2.801.gif"
        # EZC151280555F2.815.gif  EZC151280555F2.820.gif  EZC151280555F2.845.gif  EZC151280555F2.850.gif

        # specify dBZ threshold for ECHO TOP
        if values["product"] == 'EZC':
            top=chn_name.replace("'", "").replace('"', '')[-2:]
            filename=filename[:-6]+top+filename[-4:]

        print "... search for file: ", filename
        filenames=glob.glob(str(filename))
        if len(filenames) == 0:
            print "*** Error, no file found"
            quit()
        elif len(filenames) > 1:
            print "*** Warning, more than 1 datafile found: ", filenames 
        filename = filenames[0]

        print("... read data from %s" % str(filename))								 
        
        # determine the format, read data accordingly 
        from os.path import splitext
        file_name,extension = splitext(filename)

        if extension=='.gif':
            im = Image.open(str(filename))
            #im.show()

            #print 'Format (swissradar.py): ' + im.format
            #print 'Size (swissradar.py): ' + str(im.size)
            #print 'Mode (swissradar.py): ' + str(im.mode)

            radar_data = convertToValue(im, scale)
            #imo = Image.fromarray(radar_data.clip(0,255).astype("uint8"))
            #imo.save('test.png')

            print 'gif-reader: swissradar.py (min/max):', radar_data.min(), radar_data.max()

            radar_data_masked = np.ma.asarray(radar_data)
            # RZC
            #radar_data_masked.mask = (radar_data_masked == 9999.9) | (radar_data_masked <= 0.0001) 
            # VIL
            #radar_data_masked.mask = (radar_data_masked >= 9999.0) | (radar_data_masked <= 0.5) 
            # CZC
            radar_data_masked.mask = (radar_data_masked >= 9999.0) | (radar_data_masked <= 0.5)
        else:
            import sys
            sys.path.insert(0, '/opt/users/hau/PyTroll/packages/mpop/mpop/satin/metranet')
            import metranet
            ret = metranet.read_file(str(filename), physic_value=True)
            print 'metranet: swissradar.py (min/max):', ret.data.min(), ret.data.max()
            
            from numpy.ma import masked_invalid, masked_where
            radar_data_masked=masked_invalid(ret.data)
            radar_data_masked=masked_where(radar_data_masked==0,radar_data_masked)
    
        print 'swissradar.py (min/max):', radar_data_masked.min(), radar_data_masked.max()

        # save data in satscene class according to channel name
        satscene[chn_name] = radar_data_masked

        # save 3 letter radar product name 
        satscene[chn_name].product_name = values["product"]

        # save projection / area 
        #satscene.area = pyresample.utils.load_area(os.path.join(CONFIG_PATH, "areas.def"), projectionName)
        satscene[chn_name].area = pyresample.utils.load_area(os.path.join(CONFIG_PATH, "areas.def"), projectionName)

        # save units
        satscene[chn_name].units = units

def convertToValue(im, scale):
    rgb_im = im.convert('RGB')
    #r, g, b = rgb_im.getpixel((1, 1))
    #print r, g, b
   
    colorscale = np.loadtxt(scale, skiprows=1)
    #print "Colorscale: "
    #print colorscale
    translate = dict(zip(zip(colorscale[:, 1], colorscale[:,2], colorscale[:,3]), colorscale[:,-1]))

    #translate.get(rgb_im.getdata()[100], 500)  # 500 is the null value, in case the value is not found in the dictionary
    rain = np.zeros(len(rgb_im.getdata()))
    n = 0
    for i in (rgb_im.getdata()):
        #if sum(i) > 0:
        #    print i, n, translate.get(i, 9999.9)
        rain[n] = translate.get(i, 9999.9)
        n = n + 1
        

    width, height = im.size
    rain = rain.reshape(height,width)
    x = ma.masked_array(rain)
    return x

