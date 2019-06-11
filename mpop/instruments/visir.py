#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2017.

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Lars Ørum Rasmussen <ras@dmi.dk>

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

"""This module defines the generic VISIR instrument class.
"""
from mpop.imageo import geo_image
from mpop.compositer import Compositer

# pylint: disable=W0612
# remove warnings for unused prerequisites


class VisirCompositer(Compositer):

    """Compositer for Visual-IR instruments
    """

    def __call__(self, *channels, **keys):
        """Build a geoimage.
        e.g.:
        img = l.image(0.6, 0.8, -10.8, mode="RGB")
        """

        data = []
        area = None
        inv = []
        new_channels = []

        for channel in channels:
            if isinstance(channel, str):
                if channel.startswith("-"):
                    inv.append(True)
                    channel = channel[1:]
                else:
                    inv.append(False)
            else:
                if channel < 0:
                    inv.append(True)
                    channel = -channel
                else:
                    inv.append(False)

            new_channels.append(channel)

            data.append(self[channel].data)

            new_area = self[channel].area
            if area and (new_area != area):
                raise ValueError("Channels should have the same area")
            else:
                area = new_area

        self.check_channels(*new_channels)

        img = geo_image.GeoImage(data,
                                 area=area,
                                 time_slot=self.time_slot,
                                 fill_value=keys.get("fill_value", None),
                                 crange=keys.get("crange", None),
                                 mode=keys.get("mode", None))

        img.enhance(inverse=inv,
                    gamma=keys.get("gamma", 1.0),
                    stretch=keys.get("stretch", "no"))

        return img

    def channel_image(self, channel, fill_value=0):
        """Make a black and white image of the *channel*.

        Linear stretch without clipping is applied by default.
        """
        self.check_channels(channel)

        img = geo_image.GeoImage(self[channel].data,
                                 self[channel].area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="L")
        img.enhance(stretch="crude")
        return img

    def overview(self, stretch='crude', gamma=1.6, fill_value=(0, 0, 0)):
        """Make an overview RGB image composite.

        +--------------------+--------------------+
        | Channels           | Gamma (default)    |
        +====================+====================+
        | VIS0.6             | gamma 1.6          |
        +--------------------+--------------------+
        | VIS0.8             | gamma 1.6          |
        +--------------------+--------------------+
        | IR10.8 (inverted)  | gamma 1.6          |
        +--------------------+--------------------+

        Linear stretch without clipping is applied.
        """
        self.check_channels('VIS006', 'VIS008', 'IR_108')

        ch1 = self['VIS006'].check_range()
        ch2 = self['VIS008'].check_range()
        ch3 = -self['IR_108'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="RGB")

        if stretch:
            img.enhance(stretch=stretch)
        if gamma:
            img.enhance(gamma=gamma)

        return img

    overview.prerequisites = set(['VIS006', 'VIS008', 'IR_108'])

    # def overview_sun(self, stretch='crude', gamma=1.6):
    def overview_sun(self, stretch='linear', gamma=1.6, fill_value=(0, 0, 0)):
        """Make an overview RGB image composite normalising with cosine to the
        sun zenith angle.
        """
        self.check_channels('VIS006', 'VIS008', 'IR_108')

        lonlats = self['IR_108'].area.get_lonlats()

        red = self['VIS006'].sunzen_corr(self.time_slot, lonlats, limit=88.,
                                      sunmask=95).data
        green = self['VIS008'].sunzen_corr(self.time_slot, lonlats, limit=88.,
                                       sunmask=95).data
        blue = -self['IR_108'].data

        img = geo_image.GeoImage((red, green, blue),
                                 self.area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="RGB")

        if stretch:
            img.enhance(stretch=stretch)
        if gamma:
            img.enhance(gamma=gamma)

        return img

    overview_sun.prerequisites = set(['VIS006', 'VIS008', 'IR_108'])

    def night_overview(self, stretch='histogram', gamma=None):
        """Make an overview RGB image composite using IR channels.

        +--------------------+--------------------+
        | Channels           | Gamma              |
        +====================+====================+
        | IR3.9 (inverted)   | gamma 1            |
        +--------------------+--------------------+
        | IR10.8 (inverted)  | gamma 1            |
        +--------------------+--------------------+
        | IR12.0 (inverted)  | gamma 1            |
        +--------------------+--------------------+

        Histogram equalization is applied for each channel.
        """
        return self.cloudtop(stretch, gamma)

    night_overview.prerequisites = set(['IR_039', 'IR_108', 'IR_120'])

    def natural(self, stretch=None, gamma=1.8, fill_value=(0, 0, 0)):
        """Make a Natural Colors RGB image composite.

        +--------------------+--------------------+--------------------+
        | Channels           | Range (reflectance)| Gamma (default)    |
        +====================+====================+====================+
        | IR1.6              | 0 - 90             | gamma 1.8          |
        +--------------------+--------------------+--------------------+
        | VIS0.8             | 0 - 90             | gamma 1.8          |
        +--------------------+--------------------+--------------------+
        | VIS0.6             | 0 - 90             | gamma 1.8          |
        +--------------------+--------------------+--------------------+
        """
        self.check_channels('VIS006', 'VIS008', 'IR_016')

        ch1 = self['IR_016'].check_range()
        ch2 = self['VIS008'].check_range()
        ch3 = self['VIS006'].check_range()

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="RGB",
                                 crange=((0, 90),
                                         (0, 90),
                                         (0, 90)))

        if stretch:
            img.enhance(stretch=stretch)
        if gamma:
            img.enhance(gamma=gamma)

        return img

    natural.prerequisites = set(['VIS006', 'VIS008', 'IR_016'])

    def airmass(self, fill_value=(0, 0, 0)):
        """Make an airmass RGB image composite.

        +--------------------+--------------------+--------------------+
        | Channels           | Temp               | Gamma              |
        +====================+====================+====================+
        | WV6.2 - WV7.3      |     -25 to 0 K     | gamma 1            |
        +--------------------+--------------------+--------------------+
        | IR9.7 - IR10.8     |     -40 to 5 K     | gamma 1            |
        +--------------------+--------------------+--------------------+
        | WV6.2              |   243 to 208 K     | gamma 1            |
        +--------------------+--------------------+--------------------+
        """
        self.check_channels('WV_062', 'WV_073', 'IR_097', 'IR_108')

        ch1 = self['WV_062'].data - self['WV_073'].data
        ch2 = self['IR_097'].data - self['IR_108'].data
        ch3 = self['WV_062'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="RGB",
                                 crange=((-25, 0),
                                         (-40, 5),
                                         (243, 208)))
        return img

    airmass.prerequisites = set(['WV_062', 'WV_073', 'IR_097', 'IR_108'])

    def vis06(self):
        """Make a black and white image of the VIS 0.635um channel.

        Linear stretch without clipping is applied.
        """
        return self.channel_image(0.6)

    vis06.prerequisites = set(['VIS006'])

    def ir108(self):
        """Make a black and white image of the IR 10.8um channel.

        Channel is inverted. Temperature range from -70 °C (white) to
        +57.5 °C (black) is shown.
        """
        self.check_channels('IR_108')

        img = geo_image.GeoImage(self['IR_108'].data,
                                 self.area,
                                 self.time_slot,
                                 fill_value=0,
                                 mode="L",
                                 crange=(-70 + 273.15, 57.5 + 273.15))
        img.enhance(inverse=True)
        return img

    ir108.prerequisites = set(['IR_108'])

    def wv_high(self):
        """Make a black and white image of the IR 6.7um channel.

        Channel inverted and a linear stretch is applied with 0.5 %
        clipping at both ends.
        """
        self.check_channels('WV_062')

        img = geo_image.GeoImage(self['WV_062'].data,
                                 self.area,
                                 self.time_slot,
                                 fill_value=0,
                                 mode="L")
        img.enhance(inverse=True, stretch="linear")
        return img

    wv_high.prerequisites = set(['WV_062'])

    def wv_low(self):
        """Make a black and white image of the IR 7.3um channel.

        Channel data inverted and a linear stretch is applied with 0.5
        % clipping at both ends.
        """
        self.check_channels('WV_073')

        img = geo_image.GeoImage(self['WV_073'].data,
                                 self.area,
                                 self.time_slot,
                                 fill_value=0,
                                 mode="L")
        img.enhance(inverse=True, stretch="linear")
        return img

    wv_low.prerequisites = set(['WV_073'])

    def green_snow(self, fill_value=(0, 0, 0)):
        """Make a Green Snow RGB image composite.

        +--------------------+--------------------+
        | Channels           | Gamma              |
        +====================+====================+
        | IR1.6              | gamma 1.6          |
        +--------------------+--------------------+
        | VIS0.6             | gamma 1.6          |
        +--------------------+--------------------+
        | IR10.8 (inverted)  | gamma 1.6          |
        +--------------------+--------------------+

        Linear stretch without clipping.
        """
        self.check_channels('VIS006', 'IR_016', 'IR_108')

        ch1 = self['IR_016'].check_range()
        ch2 = self['VIS006'].check_range()
        ch3 = -self['IR_108'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="RGB")

        img.enhance(stretch="crude")
        img.enhance(gamma=1.6)

        return img

    green_snow.prerequisites = set(['VIS006', 'IR_016', 'IR_108'])

    def red_snow(self, fill_value=(0, 0, 0)):
        """Make a Red Snow RGB image composite.

        +--------------------+--------------------+
        | Channels           | Gamma              |
        +====================+====================+
        | VIS0.6             | gamma 1.6          |
        +--------------------+--------------------+
        | IR1.6              | gamma 1.6          |
        +--------------------+--------------------+
        | IR10.8 (inverted)  | gamma 1.6          |
        +--------------------+--------------------+

        Linear stretch without clipping.
        """
        self.check_channels('VIS006', 'IR_016', 'IR_108')

        ch1 = self['VIS006'].check_range()
        ch2 = self['IR_016'].check_range()
        ch3 = -self['IR_108'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="RGB")

        img.enhance(stretch="crude")

        return img

    red_snow.prerequisites = set(['VIS006', 'IR_016', 'IR_108'])

    def convection(self, fill_value=(0, 0, 0)):
        """Make a Severe Convection RGB image composite.

        +--------------------+--------------------+--------------------+
        | Channels           | Span               | Gamma              |
        +====================+====================+====================+
        | WV6.2 - WV7.3      |     -30 to 0 K     | gamma 1            |
        +--------------------+--------------------+--------------------+
        | IR3.9 - IR10.8     |      0 to 55 K     | gamma 1            |
        +--------------------+--------------------+--------------------+
        | IR1.6 - VIS0.6     |    -70 to 20 %     | gamma 1            |
        +--------------------+--------------------+--------------------+
        """
        self.check_channels('VIS006', 'IR_016', 'IR_039', 'WV_062', 'WV_073', 'IR_108')

        ch1 = self['WV_062'].data - self['WV_073'].data
        ch2 = self['IR_039'].data - self['IR_108'].data
        ch3 = self['IR_016'].check_range() - self['VIS006'].check_range()

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="RGB",
                                 crange=((-30, 0),
                                         (0, 55),
                                         (-70, 20)))

        return img

    convection.prerequisites = set(['VIS006', 'IR_016', 'IR_039', 'WV_062', 'WV_073', 'IR_108'])

    def dust(self, fill_value=(0, 0, 0)):
        """Make a Dust RGB image composite.

        +--------------------+--------------------+--------------------+
        | Channels           | Temp               | Gamma              |
        +====================+====================+====================+
        | IR12.0 - IR10.8    |     -4 to 2 K      | gamma 1            |
        +--------------------+--------------------+--------------------+
        | IR10.8 - IR8.7     |     0 to 15 K      | gamma 2.5          |
        +--------------------+--------------------+--------------------+
        | IR10.8             |   261 to 289 K     | gamma 1            |
        +--------------------+--------------------+--------------------+
        """
        self.check_channels('IR_087', 'IR_108', 'IR_120')

        ch1 = self['IR_120'].data - self['IR_108'].data
        ch2 = self['IR_108'].data - self['IR_087'].data
        ch3 = self['IR_108'].data
        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="RGB",
                                 crange=((-4, 2),
                                         (0, 15),
                                         (261, 289)))

        img.enhance(gamma=(1.0, 2.5, 1.0))

        return img

    dust.prerequisites = set(['IR_087', 'IR_108', 'IR_120'])

    def ash(self, fill_value=(0, 0, 0)):
        """Make a Ash RGB image composite.

        +--------------------+--------------------+--------------------+
        | Channels           | Temp               | Gamma              |
        +====================+====================+====================+
        | IR12.0 - IR10.8    |     -4 to 2 K      | gamma 1            |
        +--------------------+--------------------+--------------------+
        | IR10.8 - IR8.7     |     -4 to 5 K      | gamma 1            |
        +--------------------+--------------------+--------------------+
        | IR10.8             |   243 to 303 K     | gamma 1            |
        +--------------------+--------------------+--------------------+
        """
        self.check_channels('IR_087', 'IR_108', 'IR_120')

        ch1 = self['IR_120'].data - self['IR_108'].data
        ch2 = self['IR_108'].data - self['IR_087'].data
        ch3 = self['IR_108'].data
        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="RGB",
                                 crange=((-4, 2),
                                         (-4, 5),
                                         (243, 303)))

        return img

    ash.prerequisites = set(['IR_087', 'IR_108', 'IR_120'])

    def fog(self, fill_value=(0, 0, 0)):
        """Make a Fog RGB image composite.

        +--------------------+--------------------+--------------------+
        | Channels           | Temp               | Gamma              |
        +====================+====================+====================+
        | IR12.0 - IR10.8    |     -4 to 2 K      | gamma 1            |
        +--------------------+--------------------+--------------------+
        | IR10.8 - IR8.7     |      0 to 6 K      | gamma 2.0          |
        +--------------------+--------------------+--------------------+
        | IR10.8             |   243 to 283 K     | gamma 1            |
        +--------------------+--------------------+--------------------+
        """
        self.check_channels('IR_087', 'IR_108', 'IR_120')

        ch1 = self['IR_120'].data - self['IR_108'].data
        ch2 = self['IR_108'].data - self['IR_087'].data
        ch3 = self['IR_108'].data
        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="RGB",
                                 crange=((-4, 2),
                                         (0, 6),
                                         (243, 283)))

        img.enhance(gamma=(1.0, 2.0, 1.0))

        return img

    fog.prerequisites = set(['IR_087', 'IR_108', 'IR_120'])

    def night_fog(self, fill_value=(0, 0, 0)):
        """Make a Night Fog RGB image composite.

        +--------------------+--------------------+--------------------+
        | Channels           | Temp               | Gamma              |
        +====================+====================+====================+
        | IR12.0 - IR10.8    |     -4 to 2 K      | gamma 1            |
        +--------------------+--------------------+--------------------+
        | IR10.8 - IR3.9     |      0 to 6 K      | gamma 2.0          |
        +--------------------+--------------------+--------------------+
        | IR10.8             |   243 to 293 K     | gamma 1            |
        +--------------------+--------------------+--------------------+
        """
        self.check_channels('IR_039', 'IR_108', 'IR_120')

        ch1 = self['IR_120'].data - self['IR_108'].data
        ch2 = self['IR_108'].data - self['IR_039'].data
        ch3 = self['IR_108'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="RGB",
                                 crange=((-4, 2),
                                         (0, 6),
                                         (243, 293)))

        img.enhance(gamma=(1.0, 2.0, 1.0))

        return img

    night_fog.prerequisites = set(['IR_039', 'IR_108', 'IR_120'])

    def cloudtop(self, stretch=(0.005, 0.005), gamma=None, fill_value=(0, 0, 0)):
        """Make a Cloudtop RGB image composite.

        +--------------------+--------------------+
        | Channels           | Gamma              |
        +====================+====================+
        | IR3.9 (inverted)   | gamma 1            |
        +--------------------+--------------------+
        | IR10.8 (inverted)  | gamma 1            |
        +--------------------+--------------------+
        | IR12.0 (inverted)  | gamma 1            |
        +--------------------+--------------------+

        Linear stretch with 0.5 % clipping at both ends.
        """
        self.check_channels('IR_039', 'IR_108', 'IR_120')

        ch1 = -self['IR_039'].data
        ch2 = -self['IR_108'].data
        ch3 = -self['IR_120'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="RGB")

        if stretch:
            img.enhance(stretch=stretch)
        if gamma:
            img.enhance(gamma=gamma)

        return img

    cloudtop.prerequisites = set(['IR_039', 'IR_108', 'IR_120'])

# pylint: enable=W0612
