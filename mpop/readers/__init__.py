#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2015.

# Author(s):

#   David Hoese <david.hoese@ssec.wisc.edu>

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

"""Shared objects of the various reader classes.

"""

from mpop.plugin_base import Plugin
import logging
import numbers

LOG = logging.getLogger(__name__)


class Reader(Plugin):
    """Reader plugins. They should have a *pformat* attribute, and implement
    the *load* method. This is an abstract class to be inherited.
    """
    def __init__(self, name=None,
                 file_patterns=None,
                 filenames=None,
                 description="",
                 start_time=None,
                 end_time=None,
                 area=None,
                 sensor=None,
                 **kwargs):
        """The reader plugin takes as input a satellite scene to fill in.

        Arguments:
        - `scene`: the scene to fill.
        """
        # Hold information about channels
        self.channels = {}

        # Load the config
        Plugin.__init__(self, **kwargs)

        # Use options from the config file if they weren't passed as arguments
        self.name = self.config_options.get("name", None) if name is None else name
        self.file_patterns = self.config_options.get("file_patterns", None) if file_patterns is None else file_patterns
        self.filenames = self.config_options.get("filenames", []) if filenames is None else filenames
        self.description = self.config_options.get("description", None) if description is None else description
        self.sensor = self.config_options.get("sensor", "").split(",") if sensor is None else set(sensor)

        # These can't be provided by a configuration file
        self.start_time = start_time
        self.end_time = end_time
        self.area = area

        if self.name is None:
            raise ValueError("Reader 'name' not provided")

    def add_filenames(self, *filenames):
        self.filenames |= set(filenames)

    @property
    def channel_names(self):
        """Names of all channels configured for this reader.
        """
        return sorted(self.channels.keys())

    @property
    def sensor_names(self):
        """Sensors supported by this reader.
        """
        sensors = set()
        for chn_info in self.channels.values():
            if "sensor" in chn_info:
                sensors |= set(chn_info["sensor"].split(","))
        return sensors | self.sensor

    def load_section_reader(self, section_name, section_options):
        self.config_options = section_options

    def load_section_channel(self, section_name, section_options):
        name = section_options.get("name", section_name.split(":")[-1])
        section_options["name"] = name

        # Allow subclasses to make up their own rules about channels, but this is a good starting point
        if "file_patterns" in section_options:
            section_options["file_patterns"] = section_options["file_patterns"].split(",")
        if "wavelength_range" in section_options:
            section_options["wavelength_range"] = [float(wl) for wl in section_options["wavelength_range"].split(",")]

        self.channels[name] = section_options

    def get_channel(self, key):
        """Get the channel corresponding to *key*, either by name or centerwavelength.
        """
        # get by wavelength
        if isinstance(key, numbers.Number):
            channels = [chn for chn in self.channels.values()
                        if("wavelength_range" in chn and
                           chn["wavelength_range"][0] <= key <=chn["wavelength_range"][2])]
            channels = sorted(channels,
                              lambda ch1, ch2:
                              cmp(abs(ch1["wavelength_range"][1] - key),
                                  abs(ch2["wavelength_range"][1] - key)))

            if not channels:
                raise KeyError("Can't find any projectable at %gum" % key)
            return channels[0]
        # get by name
        else:
            return self.channels[key]

    def load(self, channels_to_load):
        """Loads the *channels_to_load* into the scene object.
        """
        raise NotImplementedError

