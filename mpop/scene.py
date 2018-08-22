#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2010-2017.

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Esben S. Nielsen <esn@dmi.dk>

# This file is part of mpop.

# mpop is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# mpop is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# mpop.  If not, see <http://www.gnu.org/licenses/>.

"""The :mod:`mpop.scene` module defines satellite scenes. They are defined as
generic classes, to be inherited when needed.

A scene is a set of :mod:`mpop.channel` objects for a given time, and sometimes
also for a given area.
"""
import ConfigParser
import copy
import datetime
import imp
import logging
import os.path
import sys
import types
import weakref

import numpy as np

import mpop.satin
from mpop import CONFIG_PATH
from mpop.channel import Channel, NotLoadedError
from mpop.utils import OrderedConfigParser

LOG = logging.getLogger(__name__)
LOG.setLevel(40)
#CRITICAL 50 #ERROR 40 #WARNING 30 #INFO 20 #DEBUG 10 #NOTSET 0

try:
    # Work around for on demand import of pyresample. pyresample depends
    # on scipy.spatial which memory leaks on multiple imports
    is_pyresample_loaded = False
    from pyresample.geometry import AreaDefinition, SwathDefinition
    import mpop.projector
    is_pyresample_loaded = True
except ImportError:
    LOG.warning("pyresample missing. Can only work in satellite projection")


class Satellite(object):

    """This is the satellite class. It contains information on the satellite.
    """

    def __init__(self, (satname, number, variant)=(None, None, None)):
        try:
            self.satname = satname or "" or self.satname
        except AttributeError:
            self.satname = satname or ""
        try:
            self.number = number or "" or self.number
        except AttributeError:
            self.number = number or ""
        try:
            self.variant = variant or "" or self.variant
        except AttributeError:
            self.variant = variant or ""

    @property
    def fullname(self):
        """Full name of the satellite, that is platform name and number
        (eg "metop02").
        """
        return self.variant + self.satname + self.number

    def sat_nr(self, string=False):
        import re
        sat_nr = re.findall(r'\d+', self.fullname)
        if len(sat_nr)>0:
          sat_nr = sat_nr[0]
        else:
          sat_nr=0
        if string:
            return sat_nr
        else:
            return int(sat_nr)

    @classmethod
    def remove_attribute(cls, name):
        """Remove an attribute from the class.
        """
        return delattr(cls, name)

    @classmethod
    def add_method(cls, func):
        """Add a method to the class.
        """
        return setattr(cls, func.__name__, func)

    def add_method_to_instance(self, func):
        """Add a method to the instance.
        """
        return setattr(self, func.__name__,
                       types.MethodType(func, self.__class__))


class SatelliteScene(Satellite):

    """This is the satellite scene class. It is a capture of the satellite
    (channels) data at given *time_slot* and *area_id*/*area*.
    """

    def __init__(self, time_slot=None, area_id=None, area=None,
                 orbit=None, satellite=(None, None, None)):

        Satellite.__init__(self, satellite)

        # if(time_slot is not None and
        #   not isinstance(time_slot, datetime.datetime)):
        #    raise TypeError("Time_slot must be a datetime.datetime instance.")

        self.time_slot = time_slot

        self.area_id = None
        self.area_def = None

        if area_id is not None:
            from warnings import warn
            warn("The *area_id* attribute is deprecated."
                 "Please use *area* instead.",
                 DeprecationWarning)
            if not isinstance(area_id, str):
                raise TypeError("Area must be a string.")

        self.area = area_id

        if area is not None:
            self.area = area

        self.orbit = orbit

        self.info = {}

        self.lat = None
        self.lon = None

    def get_area(self):
        """Getter for area.
        """
        return self.area_def or self.area_id

    def set_area(self, area):
        """Setter for area.
        """
        if area is None:
            self.area_def = None
            self.area_id = None
        elif isinstance(area, str):
            self.area_id = area
            self.area_def = None
        else:
            try:
                dummy = area.area_extent
                dummy = area.x_size
                dummy = area.y_size
                dummy = area.proj_id
                dummy = area.proj_dict
                self.area_def = area
                self.area_id = None
            except AttributeError:
                try:
                    dummy = area.lons
                    dummy = area.lats
                    self.area_def = area
                    self.area_id = None
                except AttributeError:
                    raise TypeError(("Malformed area argument. "
                                     "Should be a string or an area object. "
                                     "Not %s") % type(area))

    area = property(get_area, set_area)


class SatelliteInstrumentScene(SatelliteScene):

    """This is the satellite instrument class. It is an abstract channel
    container, from which all concrete satellite scenes should be derived.

    The constructor accepts as optional arguments the *time_slot* of the scene,
    the *area* on which the scene is defined (this can be used for slicing of
    big datasets, or can be set automatically when loading), and *orbit* which
    is a string giving the orbit number.
    """
    channel_list = []

    def __init__(self, time_slot=None, area_id=None, area=None,
                 orbit=None, satellite=(None, None, None), instrument=None):

        SatelliteScene.__init__(self, time_slot, area_id, area,
                                orbit, satellite)

        try:
            self.instrument_name = instrument or self.instrument_name
        except AttributeError:
            self.instrument_name = None

        self.channels = []

        self.end_time = None
        if isinstance(self.time_slot, (tuple, list)):
            self.time_slot, self.end_time = self.time_slot

        try:
            conf = OrderedConfigParser()
            conf.read(os.path.join(CONFIG_PATH, self.fullname + ".cfg"))

            for section in conf.sections():
                if(not section[:-1].endswith("level") and
                   not section.endswith("granules") and
                   section.startswith(self.instrument_name)):
                    name = eval(conf.get(section, "name"))
                    try:
                        w_range = eval(conf.get(section, "frequency"))
                    except ConfigParser.NoOptionError:
                        w_range = (-np.inf, -np.inf, -np.inf)
                    try:
                        resolution = eval(conf.get(section, "resolution"))
                    except ConfigParser.NoOptionError:
                        resolution = 0
                    self.channels.append(Channel(name=name,
                                                 wavelength_range=w_range,
                                                 resolution=resolution))

        except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            for name, w_range, resolution in self.channel_list:
                self.channels.append(Channel(name=name,
                                             wavelength_range=w_range,
                                             resolution=resolution))

        self.channels_to_load = set([])

    def __getitem__(self, key, aslist=False):
        if isinstance(key, float):
            channels = [chn for chn in self.channels
                        if(hasattr(chn, "wavelength_range") and
                           chn.wavelength_range[0] <= key and
                           chn.wavelength_range[2] >= key)]
            channels = sorted(channels,
                              lambda ch1, ch2:
                                  ch1.__cmp__(ch2, key))

        elif isinstance(key, str):
            channels = [chn for chn in self.channels
                        if chn.name == key]
            channels = sorted(channels)

        elif isinstance(key, int):
            channels = [chn for chn in self.channels
                        if int(np.round(chn.resolution)) == key]
            channels = sorted(channels)

        elif isinstance(key, (tuple, list)):
            if len(key) == 0:
                raise KeyError("Key list must contain at least one element.")
            channels = self.__getitem__(key[0], aslist=True)
            if len(key) > 1 and len(channels) > 0:
                dummy_instance = SatelliteInstrumentScene()
                dummy_instance.channels = channels
                channels = dummy_instance.__getitem__(key[1:], aslist=True)
        else:
            raise TypeError("Malformed key: " + str(key))

        if len(channels) == 0:
            raise KeyError("No channel corresponding to " + str(key) + ".")
        elif aslist:
            return channels
        else:
            return channels[0]

    def __delitem__(self, key):
        """Delete a channel.
        """
        chan = self[key]
        self.channels.remove(chan)
        del chan

    def __setitem__(self, key, data):
        # Add a channel if it is not already in the scene. Works only if key is
        # a string.
        try:
            if key not in self:
                # if it's a blob with name and data, add it as is.
                if hasattr(data, "name") and hasattr(data, "data"):
                    self.channels.append(data)
                else:
                    kwargs = {"name": key}
                    for attr in ["wavelength_range", "resolution"]:
                        try:
                            kwargs[attr] = getattr(data, attr)
                        except (AttributeError, NameError):
                            pass
                    self.channels.append(Channel(**kwargs))
        except AttributeError:
            pass

        # Add the data.
        if isinstance(data, np.ma.core.MaskedArray):
            self[key].data = data
        else:
            try:
                self[key].data = data.data
            except AttributeError:
                self[key].data = data

        # if isinstance(data, Channel):
        #     self.channels.append(Channel(name=key,
        #                              wavelength_range=data.wavelength_range,
        #                              resolution=data.resolution))
        #     self[key].data = data.data
        # else:
        #     try:
        #         self[key].data = data
        #     except KeyError:
        #         self.channels.append(Channel(name=key))
        #         self[key].data = data

    def __str__(self):
        return "\n".join([str(chn) for chn in self.channels])

    def __iter__(self):
        return self.channels.__iter__()

    def _set_reader(self, pformat):
        """Gets the reader for *pformat* format, and puts it in the `reader`
        attribute.
        """

        elements = pformat.split(".")
        if len(elements) == 1:
            reader_module = pformat
            reading_element = 'load'

            reader = "mpop.satin." + reader_module

            # Loading old style plugins
            reader_module = pformat
            LOG.info("old style plugin: " + pformat)

            try:
                # Look for builtin reader
                imp.find_module(reader_module, mpop.satin.__path__)
            except ImportError:
                # Look for custom reader
                loader = __import__(reader_module, globals(),
                                    locals(), [reading_element])
            else:
                loader = __import__(reader, globals(),
                                    locals(), [reading_element])

            # Build a custom Reader plugin on the fly...
            from mpop.plugin_base import Reader
            reader_class = type(elements[-1].capitalize() + "Reader",
                                (Reader,),
                                {"pformat": elements[-1]})

            reader_instance = reader_class(self)

            # ... and set its "load" attribute with the "load" function of the
            # loader module
            loader = getattr(loader, reading_element)
            setattr(reader_instance, reading_element, loader)

            setattr(self, elements[-1] + "_reader", reader_instance)

        else:
            reader_module = ".".join(elements[:-1])
            reading_element = elements[-1]

            reader = "mpop.satin." + reader_module
            try:
                # Look for builtin reader
                imp.find_module(reader_module, mpop.satin.__path__)
            except ImportError:
                # Look for custom reader
                loader = __import__(reader_module, globals(),
                                    locals(), [reading_element])
            else:
                loader = __import__(reader, globals(),
                                    locals(), [reading_element])

            loader = getattr(loader, reading_element)
            reader_instance = loader(self)
            setattr(self, loader.pformat + "_reader", reader_instance)

        return reader_instance

    def load(self, channels=None, load_again=False, area_extent=None, **kwargs):
        """Load instrument data into the *channels*. *Channels* is a list or a
        tuple containing channels we will load data into, designated by there
        center wavelength (float), resolution (integer) or name (string). If
        None, all channels are loaded.

        The *load_again* boolean flag allows to reload the channels even they
        have already been loaded, to mirror changes on disk for example. This
        is false by default.

        The *area_extent* keyword lets you specify which part of the data to
        load. Given as a 4-element sequence, it defines the area extent to load
        in satellite projection.

        The other keyword arguments are passed as is to the reader
        plugin. Check the corresponding documentation for more details.
        """

        # Set up the list of channels to load.
        if channels is None:
            for chn in self.channels:
                self.channels_to_load |= set([chn.name])

        elif isinstance(channels, (list, tuple, set)):
            self.channels_to_load = set()
            for chn in channels:
                try:
                    self.channels_to_load |= set([self[chn].name])
                except KeyError:
                    self.channels_to_load |= set([chn])

        else:
            raise TypeError("Channels must be a list/"
                            "tuple/set of channel keys!")

        loaded_channels = [chn.name for chn in self.loaded_channels()]
        if load_again:
            for chn in self.channels_to_load:
                if chn in loaded_channels:
                    self.unload(chn)
                    loaded_channels = []
        else:
            for chn in loaded_channels:
                self.channels_to_load -= set([chn])
                
        # find the plugin to use from the config file
        conf = ConfigParser.ConfigParser()
        try:
            conf.read(os.path.join(CONFIG_PATH, self.fullname + ".cfg"))
            if len(conf.sections()) == 0:
                raise ConfigParser.NoSectionError(("Config file did "
                                                   "not make sense"))
            levels = [section for section in conf.sections()
                      if section.startswith(self.instrument_name + "-level")]
        except ConfigParser.NoSectionError:
            LOG.warning("Can't load data, no config file for " + self.fullname)
            self.channels_to_load = set()
            return

        levels.sort()

        # skip level1 -> reserved for compressed HRIT
        if levels[0] == self.instrument_name + "-level1":
            levels = levels[1:]

        if len(levels) == 0:
            raise ConfigParser.NoSectionError(
                self.instrument_name + "-levelN (N>1) to tell me how to" +
                " read data... Not reading anything.")

        for level in levels:
            if len(self.channels_to_load) == 0:
                return

            if "reader_level" in kwargs.keys():
                if kwargs["reader_level"] is not None:
                    if kwargs["reader_level"] not in levels:
                        LOG.warning("Specified reader_level \""+
                                    kwargs["reader_level"]+"\" in not defined in config file")
                        print "*** ERROR, specified reader_level \""+\
                                    kwargs["reader_level"]+"\" in not defined in config file"
                    LOG.debug(
                        "Using explecit definition of reader level: " + kwargs["reader_level"])
                    if kwargs["reader_level"] != level:
                        continue

            LOG.debug("Looking for sources in section " + level)
            reader_name = conf.get(level, 'format')
            try:
                reader_name = eval(reader_name)
            except NameError:
                reader_name = str(reader_name)
            LOG.debug("Using plugin mpop.satin." + reader_name)

            # read the data
            reader = "mpop.satin." + reader_name
            try:
                reader_instance = self._set_reader(reader_name)
                if area_extent is not None:
                    if(isinstance(area_extent, (tuple, list)) and
                       len(area_extent) == 4):
                        kwargs["area_extent"] = area_extent
                    else:
                        raise ValueError("Area extent must be a sequence of "
                                         "four numbers.")

                reader_instance.load(self, **kwargs)
            except ImportError, err:
                LOG.exception("ImportError while loading " + reader_name + ": "
                              + str(err))
                continue
            loaded_channels = set([chn.name for chn in self.loaded_channels()])
            just_loaded = loaded_channels & self.channels_to_load
            if len(just_loaded) == 0:
                LOG.info("No channels loaded with " + reader + ".")
            self.channels_to_load -= loaded_channels
            LOG.debug("Successfully loaded: " + str(just_loaded))

        if len(self.channels_to_load) > 0:
            LOG.warning("Unable to import channels "
                        + str(self.channels_to_load))

        self.channels_to_load = set()

    def save(self, filename, to_format="netcdf4", **options):
        """Saves the current scene into a file of format *to_format*. Supported
        formats are:

        - *netcdf4*: NetCDF4 with CF conventions.
        """

        writer = "satout." + to_format
        try:
            writer_module = __import__(writer, globals(), locals(), ["save"])
        except ImportError, err:
            raise ImportError("Cannot load " + writer + " writer: " + str(err))

        return writer_module.save(self, filename, **options)

    def unload(self, *channels):
        """Unloads *channels* from
        memory. :meth:`mpop.scene.SatelliteInstrumentScene.load` must be called
        again to reload the data.
        """
        for chn in channels:
            try:
                self[chn].data = None
            except AttributeError:
                LOG.warning("Can't unload channel" + str(chn))

    def add_to_history(self, message):
        """Adds a message to history info.
        """
        timestr = datetime.datetime.utcnow().isoformat()
        timed_message = str(timestr + " - " + message)
        if not self.info.get("history", ""):
            self.info["history"] = timed_message
        else:
            self.info["history"] += "\n" + timed_message

    def check_channels(self, *channels):
        """Check if the *channels* are loaded, raise an error otherwise.
        """
        for chan in channels:
            if not self[chan].is_loaded():
                raise NotLoadedError("Required channel %s not loaded,"
                                     " aborting." % chan)

        return True

    def loaded_channels(self):
        """Return the set of loaded_channels.
        """
        return set([chan for chan in self.channels if chan.is_loaded()])

    def get_orbital(self, allow_NEAR_NORM=False):
        from pyorbital.orbital import Orbital
        from pyorbital import tlefile

        from pyorbital.tlefile import get_norad_line
        sat_line = get_norad_line(self.satname, self.number)
        self.orbital = Orbital(sat_line, allow_NEAR_NORM=allow_NEAR_NORM) 

        return self.orbital

    def estimate_cth(self, cth_atm="best", time_slot=None, IR_108=None):
        """
        General purpose
        ===============
           Estimation of the cloud top height using the 10.8 micron channel
           limitations: this is the most simple approach
           a simple fit of the ir108 to the temperature profile
                * no correction for water vapour or any other trace gas
                * no viewing angle dependency
                * no correction for semi-transparent clouds
                * no special treatment of temperature inversions
        Example call
        ============
           data.estimate_cth(cth_atm="best")

        input arguments
        ===============
          cth_atm    * using temperature profile to estimate the cloud top height
                       possible choices are (see estimate_cth in mpop/tools.py):
                       "standard", "tropics", "midlatitude summer", "midlatitude winter", "subarctic summer", "subarctic winter"
                       this will choose the corresponding atmospheric AFGL temperature profile
                     * new choice: "best" -> choose according to central (lon,lat) and time_slot from:
                       "tropics", "midlatitude summer", "midlatitude winter", "subarctic summer", "subarctic winter"
          time_slot  * current observation time as (datetime.datetime() object)
                       used for choosing most represenative standard temperature profile 
                       time_slot option can be omitted, the function tries to use self.time_slot
          IR_108     * explicit IR_108 channel, 
                       handy, if IR_108 is not loaded into the same satellite scene object 
                       (e.g. because it comes from another satellite)
        """

        print "*** Simple estimation of Cloud Top Height with IR_108 channel"

        # check if IR_108 is loaded

        if IR_108 is not None:
            if isinstance(IR_108, mpop.channel.Channel):
                ir108 = IR_108.data
            else:
                print "*** Error in estimate_cth (mpop.scene.py)"
                print "    given IR_108 is not of the type mpop.channel.Channel"
                quit()
        else:
            loaded_channels = [chn.name for chn in self.loaded_channels()]
            if "IR_108" in loaded_channels:
                IR_108 = self["IR_108"]
                ir108 = self["IR_108"].data
            else:
                print "*** Error in estimate_cth (mpop/scene.py)"
                print "    IR_108 is required to estimate CTH, but not loaded"
                quit()

        # choose atmosphere
        if cth_atm.lower() == "best":
            # get central lon/lat coordinates
            (yc, xc) = ir108.shape
            (lon, lat) = self.area.get_lonlat(yc / 2, xc / 2)

            if time_slot is None:
                if hasattr(self, 'time_slot'):
                    time_slot = self.time_slot
                else:
                    print "*** Error, in estimate_cth (mpop/channel.py)"
                    print "    when using cth_atm=\"best\" also the time_slot information is required!"
                    quit()

            # automatic choise of temperature profile
            doy = time_slot.timetuple().tm_yday
            print "... automatic choise of temperature profile lon=", lon, " lat=", lat, ", time=", str(time_slot), ", doy=", doy
            if abs(lat) <= 30.0:
                cth_atm = "tropics"
            elif doy < 80 or doy <= 264:
                # northern summer
                if lat < -60.0:
                    cth_atm = "subarctic winter"
                elif -60.0 <= lat and lat < -30.0:
                    cth_atm = "midlatitude winter"
                elif 30.0 < lat and lat <= 60.0:
                    cth_atm = "midlatitude summer"
                elif 60.0 < lat:
                    cth_atm = "subarctic summer"
            else:
                # northern winter
                if lat < -60.0:
                    cth_atm = "subarctic summer"
                elif -60.0 <= lat and lat < -30.0:
                    cth_atm = "midlatitude summer"
                elif 30.0 < lat and lat <= 60.0:
                    cth_atm = "midlatitude winter"
                elif 60 < lat:
                    cth_atm = "subarctic winter"
            print "    choosing temperature profile for ", cth_atm

        # estimate cloud top height by searching first fit of ir108 with
        # temperature profile
        from mpop.tools import estimate_cth
        cth = estimate_cth(ir108, cth_atm=cth_atm)

        # create new channel named "CTH"
        self.channels.append(Channel(name="CTH",
                                     wavelength_range=[0., 0., 0.],
                                     resolution=self["IR_108"].resolution,
                                     data=cth,
                                     calibration_unit="m"))

        # copy additional information from IR_108
        self["CTH"].info = IR_108.info
        self["CTH"].info['units'] = 'm'
        self["CTH"].area = IR_108.area
        self["CTH"].area_id = IR_108.area_id
        self["CTH"].area_def = IR_108.area_def
        self["CTH"].resolution = IR_108.resolution

        return cth

    def parallax_corr(self, fill="False", estimate_cth=False, cth=None, IR_108=None, cth_atm='best', time_slot=None, replace=False):
        """
        General purpose
        ===============
           perform the CTH parallax corretion for all loaded channels
 
       Example call
        ============
           global_sat = global_sat.parallax_corr(fill='bilinear', estimate_cth=True, replace=True)
        input arguments
        ===============
          fill          * string: when shifing gaps occur in the dataset
                          if "False" -> dont fill gaps
                          otherwise give method for interpolation, e.g. "nearest" or "bilinear"
                          see channel.py -> function parallax_corr 
          estimate_cth  * estimate cloud top height or use CTH from NWC-SAF
                          either assumes that CTTH is loaded into the global_data 
                          or that cth is given explicitely by the cth argument
          cth           * cloud top height masked array 
                          has to be in exactly same projection as all loaded channels
          cth_atm       * using temperature profile to estimate the cloud top height
                          see estimate_cth
          time_slot     * auxilary information to choose the most representative 
                          standard temperature profile
          replace       * boolean: if true, data inside channel object is replaced
                          if false, new channel created, e.g. "IR_108" -> "IR_108_PC"
        """ 

        loaded_channels = [chn.name for chn in self.loaded_channels()]
        if len(loaded_channels) == 0:
            return

        # currently orbit calculation for deep space (=geostationary) satellites are not yet implemented
        # we use near space calculations instead, !!! PLEASE FIX IF YOU HAVE TIME !!!
        allow_NEAR_NORM=True
        
        # loop over channels and check, if one is a normal radiance channel
        # having the method to calculate the viewing geometry
        for chn in self.loaded_channels():
            if hasattr(chn, 'get_viewing_geometry'):
                # calculate the viewing geometry of the SEVIRI sensor
                print "... calculate viewing geometry using ", chn.name
                (azi, ele) = chn.get_viewing_geometry(self.get_orbital(allow_NEAR_NORM=allow_NEAR_NORM), self.time_slot)
                break

        # choose best way to get CTH for parallax correction
        if not estimate_cth:
            if cth is None:
                if "CTTH" in loaded_channels:
                    # make a copy of CTH, as it might get replace by its parallax
                    # corrected version
                    cth = copy.deepcopy(self["CTTH"].height)
                elif "CTH" in loaded_channels:
                    # make a copy of CTH, as it might get replace by its parallax
                    # corrected version
                    cth = copy.deepcopy(self["CTH"].data)
                else:
                    print "*** Error in parallax_corr (mpop.scene.py)"
                    print "    parallax correction needs some cloud top height information"
                    print "    please load the NWC-SAF CTTH product (recommended) or"
                    print "    activate the option data.parallax_corr( estimate_cth=True )"
                    quit()
            else:
                # check if cth (given as argument to this function) has the same size as all loaded channels
                for chn in self.loaded_channels():
                    if isinstance(chn, mpop.channel.Channel):
                        #print type(chn)
                        if not chn.data.shape==cth.shape:
                            raise ValueError("The shape of the channels and the cth must agree.")
 
        else:
            if IR_108 is not None:
                # estimate CTH with explicit IR_108
                print "*** Warning, it is not recommended to use IR_108 observations from another platform"
                print "    the IR_108 (with parallax shift) might represent a neighbouring cloud, if viewing geometry differs"
                self.estimate_cth(cth_atm=cth_atm, time_slot=time_slot, IR_108=IR_108)
            else:
                if "IR_108" in loaded_channels:
                    # estimate CTH with IR_108 that is included in self
                    self.estimate_cth(cth_atm=cth_atm, time_slot=time_slot)
                else:
                    print "*** Error in parallax_corr (mpop.scene.py)"
                    print "    parallax correction needs some cloud top height information."
                    print "    You specified the estimation of CTH with the IR_108, "
                    print "    but IR_108 channel is not loaded"
                    quit()
                    
            # copy cloud top height data
            cth = self["CTH"].data
        
        # perform parallax correction for each loaded channel
        for chn in self.loaded_channels():
            if isinstance(chn, mpop.channel.Channel):
                if hasattr(chn, 'parallax_corr'):
                    print "... perform parallax correction for \'"+ chn.name+"\'" #, chn.area
                    if replace:
                        chn_name_PC = chn.name
                        print "    replace channel ", chn_name_PC
                    else:
                        chn_name_PC = chn.name + "_PC"
                        print "    create channel ", chn_name_PC

                    # take care of the parallax correction
                    self[chn_name_PC] = chn.parallax_corr(cth=cth, azi=azi, ele=ele, fill=fill)
                else:
                    LOG.warning("Channel " + str(chn.name) + " has no attribute parallax_corr,"
                                "thus parallax effect wont be corrected.")
                    print "Channel " + str(chn.name) + " has no attribute parallax_corr,"
                    print "thus parallax effect wont be corrected."
            else:
                print "... Warning, cannot perform parallax correction for tpye: ", type(chn)

        return self

    def project(self, dest_area, channels=None, precompute=False, mode=None,
                radius=None, nprocs=1):
        """Make a copy of the current snapshot projected onto the
        *dest_area*. Available areas are defined in the region configuration
        file (ACPG). *channels* tells which channels are to be projected, and
        if None, all channels are projected and copied over to the return
        snapshot.

        If *precompute* is set to true, the projecting data is saved on disk
        for reusage. *mode* sets the mode to project in: 'quick' which works
        between cartographic projections, and, as its denomination indicates,
        is quick (but lower quality), and 'nearest' which uses nearest
        neighbour for best projection. A *mode* set to None uses 'quick' when
        possible, 'nearest' otherwise.

        *radius* defines the radius of influence for neighbour search in
        'nearest' mode (in metres). Setting it to None, or omitting it will
        fallback to default values (5 times the channel resolution) or 10,000m
        if the resolution is not available.

        Note: channels have to be loaded to be projected, otherwise an
        exception is raised.
        """

        if not is_pyresample_loaded:
            # Not much point in proceeding then
            return self

        _channels = set([])

        if channels is None:
            for chn in self.loaded_channels():
                _channels |= set([chn])

        elif isinstance(channels, (list, tuple, set)):
            for chn in channels:
                try:
                    _channels |= set([self[chn]])
                except KeyError:
                    LOG.warning("Channel " + str(chn) + " not found,"
                                "thus not projected.")
        else:
            raise TypeError("Channels must be a list/"
                            "tuple/set of channel keys!")

        res = copy.copy(self)

        if isinstance(dest_area, str):
            dest_area = mpop.projector.get_area_def(dest_area)

        res.area = dest_area
        res.channels = []

        if not _channels <= self.loaded_channels():
            LOG.warning("Cannot project nonloaded channels: %s.",
                        _channels - self.loaded_channels())
            LOG.info("Will project the other channels though.")
            _channels = _channels and self.loaded_channels()

        cov = {}

        for chn in sorted(_channels, key=lambda x: x.resolution, reverse=True):
            if chn.area is None:
                if self.area is None:
                    area_name = ("swath_" + self.fullname + "_" +
                                 str(self.time_slot) + "_"
                                 + str(chn.shape))
                    chn.area = area_name
                else:
                    if is_pyresample_loaded:
                        try:
                            chn.area = AreaDefinition(
                                self.area.area_id + str(chn.shape),
                                self.area.name,
                                self.area.proj_id,
                                self.area.proj_dict,
                                chn.shape[1],
                                chn.shape[0],
                                self.area.area_extent,
                                self.area.nprocs)

                        except AttributeError:
                            try:
                                dummy = self.area.lons
                                dummy = self.area.lats
                                chn.area = self.area
                                area_name = ("swath_" + self.fullname + "_" +
                                             str(self.time_slot) + "_"
                                             + str(chn.shape))
                                chn.area.area_id = area_name
                            except AttributeError:
                                chn.area = self.area + str(chn.shape)
                    else:
                        chn.area = self.area + str(chn.shape)
            else:  # chn.area is not None
                # if (is_pyresample_loaded and
                #     (not hasattr(chn.area, "area_id") or
                #      not chn.area.area_id)):
                #     area_name = ("swath_" + self.fullname + "_" +
                #                  str(self.time_slot) + "_"
                #                  + str(chn.shape) + "_"
                #                  + str(chn.name))
                #     chn.area.area_id = area_name

                # This leaks memory !
                #LOG.debug("chn.area = " + str(chn.area))
                #LOG.debug("type(chn.area) = " + str(type(chn.area)))
                if is_pyresample_loaded:
                    area_name = ("swath_" + self.fullname + "_" +
                                 str(self.time_slot) + "_" +
                                 str(chn.shape) + "_" +
                                 str(chn.name))
                    LOG.debug("pyresample is loaded... area-name = " +
                              str(area_name))
                    if hasattr(chn.area, "area_id") and not chn.area.area_id:
                        LOG.debug("chn.area has area_id attribute...")
                        chn.area.area_id = area_name
                    elif not hasattr(chn.area, "area_id") and not isinstance(chn.area, str):
                        setattr(chn.area, 'area_id', area_name)

            if isinstance(chn.area, str):
                area_id = chn.area
            else:
                area_id = chn.area_id or chn.area.area_id

            if area_id not in cov:
                if radius is None:
                    if chn.resolution > 0:
                        radius = 5 * chn.resolution
                    else:
                        radius = 10000

                cov[area_id] = mpop.projector.Projector(chn.area,
                                                        dest_area,
                                                        mode=mode,
                                                        radius=radius,
                                                        nprocs=nprocs)
                if precompute:
                    try:
                        cov[area_id].save()
                    except IOError:
                        LOG.exception("Could not save projection.")

            try:
                res.channels.append(chn.project(cov[area_id]))
            except NotLoadedError:
                LOG.warning("Channel " + str(chn.name) + " not loaded, "
                            "thus not projected.")

        # Compose with image object
        try:
            if res._CompositerClass is not None:
                # Pass weak ref to compositor to allow garbage collection
                res.image = res._CompositerClass(weakref.proxy(res))
        except AttributeError:
            pass

        return res

if sys.version_info < (2, 5):
    def any(iterable):
        for element in iterable:
            if element:
                return True
        return False


def assemble_segments(segments):
    """Assemble the scene objects listed in *segment_list* and returns the
    resulting scene object.
    """
    from mpop.satellites import GenericFactory

    channels = set([])
    for seg in segments:
        channels |= set([chn.name for chn in seg.loaded_channels()])

    seg = segments[0]

    new_scene = GenericFactory.create_scene(seg.satname, seg.number,
                                            seg.instrument_name, seg.time_slot,
                                            seg.orbit, variant=seg.variant)

    swath_definitions = {}

    for chn in channels:
        new_scene[chn] = np.ma.concatenate([seg[chn].data
                                            for seg in segments
                                            if seg[chn].is_loaded()])
        try:

            area_names = tuple([seg[chn].area.area_id
                                for seg in segments
                                if seg[chn].is_loaded()])
            if area_names not in swath_definitions:

                lons = np.ma.concatenate([seg[chn].area.lons[:]
                                          for seg in segments
                                          if seg[chn].is_loaded()])
                lats = np.ma.concatenate([seg[chn].area.lats[:]
                                          for seg in segments
                                          if seg[chn].is_loaded()])
                new_scene[chn].area = SwathDefinition(lons=lons, lats=lats)
                area_name = "+".join(area_names)
                new_scene[chn].area.area_id = area_name
                new_scene[chn].area_id = area_name
                swath_definitions[area_names] = new_scene[chn].area
            else:
                new_scene[chn].area = swath_definitions[area_names]
                new_scene[chn].area_id = new_scene[chn].area.area_id
        except AttributeError:
            pass

    try:
        lons = np.ma.concatenate([seg.area.lons[:] for seg in segments])
        lats = np.ma.concatenate([seg.area.lats[:] for seg in segments])
        new_scene.area = SwathDefinition(lons=lons, lats=lats)
        for chn in channels:
            if any([seg[chn].area for seg in segments]):
                try:
                    lon_arrays = []
                    lat_arrays = []
                    for seg in segments:
                        if seg[chn].area is not None:
                            lon_arrays.append(seg[chn].area.lons[:])
                            lat_arrays.append(seg[chn].area.lats[:])
                        else:
                            lon_arrays.append(seg.area.lons[:])
                            lat_arrays.append(seg.area.lats[:])
                    lons = np.ma.concatenate(lon_arrays)
                    lats = np.ma.concatenate(lat_arrays)
                    new_scene[chn].area = SwathDefinition(lons=lons, lats=lats)
                except AttributeError:
                    pass
    except AttributeError:
        pass

    return new_scene
