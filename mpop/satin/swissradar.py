from __future__ import division
from __future__ import print_function

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

        for i in range(9):
            radar_product='radar-{0:1d}'.format(i+1)
            prod_name = conf.get(radar_product, "name")
            if prod_name.replace("'", "").replace('"', '') == chn_name.replace("'", "").replace('"', ''):
                values["product"]=conf.get(radar_product, "product_name").replace("'", "").replace('"', '')
                if verbose:
                    print("... radar product to read:", chn_name, ' (', values["product"], ')')
                scale = conf.get(radar_product, "scale")
                if verbose:
                    print('... read scale from: ', scale)
                units = conf.get(radar_product, "units").replace("'", "").replace('"', '')
                print(chn_name, units)
                break

        filename = os.path.join(
            satscene.time_slot.strftime(conf.get("radar-level2", "dir", raw=True)) % values,
            satscene.time_slot.strftime(conf.get("radar-level2", "filename", raw=True)) % values)

        # Load data from file, possible filenames 
        # meteoswiss.radar.precip.201906071122.gif
        # EZC151280555F2.815.gif  EZC151280555F2.820.gif  EZC151280555F2.845.gif  EZC151280555F2.850.gif
        # EZC181501700VL.815      EZC151280555F2.820      EZC151280555F2.845      EZC151280555F2.850
        # RZC143081530F2.801.gif

        if filename[:10]=="meteoswiss" and values["product"] != "RZC":
            print("*** ERROR in load (swissradar.py)")
            print("  for near real time only the product precip (RZC) is implemented")
            quit()

        # determine the format, read data accordingly 
        from os.path import splitext
        file_name,extension = splitext(filename)
            
        # specify dBZ threshold for ECHO TOP
        if values["product"] == 'EZC':
            #chn_name='EchoTOP15'
            top=chn_name.replace("'", "").replace('"', '')[-2:]
            print("ECHOTOP", top)
            #%(product)s%y%j%H%M??.???*
            if extension=='.gif':
                filename=filename[:-3]+top+filename[-1:]
            else:
                filename=filename[:-3]+top+filename[-1:]

        print("... search for file: ", filename)
        filenames=glob.glob(str(filename))
        if len(filenames) == 0:
            print("*** Error, no file found")
            quit()
        elif len(filenames) > 1:
            print("*** Warning, more than 1 datafile found: ", filenames)
            ## for echotop select the correct file
            #if values["product"] == 'EZC':
            #    print ("*** TO BE IMPLEMENTED: Choose correct file for EchoTOP", filenames)                 
            #else:
            #    print ("*** Warning, more than 1 datafile found: ", filenames)
        filename = filenames[0]

        print(("... read data from %s" % str(filename)))								 
        
        if extension=='.gif':
            im = Image.open(str(filename))
            #im.show()

            #print 'Format (swissradar.py): ' + im.format
            #print 'Size (swissradar.py): ' + str(im.size)
            #print 'Mode (swissradar.py): ' + str(im.mode)

            radar_data = convertToValue(im, scale)
            #imo = Image.fromarray(radar_data.clip(0,255).astype("uint8"))
            #imo.save('test.png')

            print('gif-reader: swissradar.py (min/max):', radar_data.min(), radar_data.max())

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
            from . import metranet
            ret = metranet.read_file(str(filename), physic_value=True)
            #print dir(ret)
            ## ['__doc__', '__init__', '__module__', 'data', 'header', 'moment', 'pol_header', 'scale', 'type']
            #print ret.__doc__
            #print ret.header
            ##{u'usr_radar_required': u'ADLPW', u'data_unit': u'kg/m2', u'pid': u'LZC', u'table_endianness': u'big',
            ## u'usr_quality': u'55555', u'total_sweep': u'20', u'usr_visib_upper': u'200', u'rect_yres': u'1.000000',
            ## u'quality': u'77777', u'table_size': u'1024', u'data_width': u'0.499998', u'volume_time': u'1527696000',
            ## u'p_ver': u'1.0.2', u'table_num': u'1', u'compressed_bytes': u'8044', u'usr_zwB': u'0.57', u'usr_zwA': u'3.44',
            ## 'row': u'640', u't_ver': u'2.4.3.24 Oct  9 2017', u'usr_visib_lower': u'50', u'product': u'VIL', u'data_type': u'BYTE',
            ## u'format': u'RECT', u'moment': u'UZ', u'uncompressed_bytes': u'454400', u'data_bits': u'8', u'usr_max_dbz': u'57.00',
            ## u'a_ver': u'1.0.0', u'usr_speckle_filt': u'0', u'radar': u'ADLPW', 'column': u'710', u'usr_min_height': u'1.50',
            ## u'rect_xres': u'1.000000', u'table_name': u'8bit_metranet_vil_0.5res', u'time': u'1815016000'}
            #print type(ret.data)
            ##<type 'numpy.ndarray'>
            print('metranet B: swissradar.py (min/max):', np.nanmin(ret.data), np.nanmax(ret.data))
            
            from numpy.ma import masked_invalid, masked_where
            radar_data_masked=masked_invalid(ret.data)
            radar_data_masked=np.nan_to_num(radar_data_masked)  # replace nan with 0, and inf with large number 
            radar_data_masked=masked_where(radar_data_masked==0,     radar_data_masked)
            radar_data_masked=masked_where(radar_data_masked==9999.0,radar_data_masked)
            ##radar_data_masked=masked_where(radar_data_masked<0.05, radar_data_masked)
            radar_data_masked.mask = (radar_data_masked >= 9999.0) | (radar_data_masked <= 0.5)
            print(type(radar_data_masked))
    
        #print 'swissradar.py (min/max):', radar_data_masked.min(), radar_data_masked.max()

        # save data in satscene class according to channel name
        satscene[chn_name] = radar_data_masked
        #satscene[chn_name] = ret.data

        # save 3 letter radar product name 
        satscene[chn_name].product_name = values["product"]

        # save projection / area 
        #satscene.area = pyresample.load_area(os.path.join(CONFIG_PATH, "areas.def"), projectionName)
        satscene[chn_name].area = pyresample.load_area(os.path.join(CONFIG_PATH, "areas.def"), projectionName)

        # copy units from metranet header to pytroll object
        if extension=='.gif':
            satscene[chn_name].units = units
        else:
            satscene[chn_name].units = ret.header['data_unit']
        
def convertToValue(im, scale):
    rgb_im = im.convert('RGB')
    #r, g, b = rgb_im.getpixel((1, 1))
    #print r, g, b
   
    colorscale = np.loadtxt(scale, skiprows=1)
    #print ("Colorscale: ")
    #print (colorscale)
    translate = dict(list(zip(list(zip(colorscale[:, 1], colorscale[:,2], colorscale[:,3])), colorscale[:,-1])))

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

