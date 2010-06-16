#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2010.

# SMHI,
# Folkborgsvägen 1,
# Norrköping, 
# Sweden

# Author(s):
 
#   Martin Raspaud <martin.raspaud@smhi.se>

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

"""This module is the loader for Metop02 scenes, ears (eumetcast) version.
"""

from pp.instruments.avhrr import AvhrrScene


class EarsMetop02AvhrrScene(AvhrrScene):
    """This class implements Metop02 scenes as captured by the Avhrr
    instrument. It's constructor accepts the same arguments as
    :class:`pp.scene.SatelliteScene`.
    """
    
    satname = "metop"
    number = "02"
    variant = "ears"
    
    lat = None
    lon = None

