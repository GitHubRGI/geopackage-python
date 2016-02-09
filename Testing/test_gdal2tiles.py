#!/usr/bin/python2.7
"""
Copyright (C) 2014 Reinventing Geospatial, Inc.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>,
or write to the Free Software Foundation, Inc., 59 Temple Place -
Suite 330, Boston, MA 02111-1307, USA.

Authors:
    Steven D. Lander
Date: 2015-11
Requires: Python Imaging Library (PIL or Pillow), Sqlite3
Description: Test cases for gdal2tiles_parallel.py

Version:
"""

from os.path import abspath
from sys import path

path.append(abspath("Tiling"))
from gdal2tiles_parallel import Tile
from gdal2tiles_parallel import GlobalMercatorProfile


class TestGlobalMercatorProfile:
    """Test the GlobalMercatorProfile object"""

    def test_tile_size(self):
        """Test the default tile size"""
        gmp = GlobalMercatorProfile()
        assert gmp.tile_size == 256

    def test_tile_size_custom(self):
        """Test a custom tile size"""
        gmp = GlobalMercatorProfile(512)
        assert gmp.tile_size == 512

    def test_lower_left_tile(self):
        """Test lower left tile method"""
        gmp = GlobalMercatorProfile()
        tile = Tile(0, 1)
        zoom = 1
        result_tile, result_zoom = gmp.lower_left_tile(tile, zoom)
        assert result_tile.tx == tile.tx and \
            result_tile.ty == tile.ty and \
            result_zoom == zoom

    def test_upper_left_tile(self):
        """Test upper left tile method"""
        gmp = GlobalMercatorProfile()
        tile = Tile(0, 1)
        zoom = 1
        result_tile, result_zoom = gmp.upper_left_tile(tile, zoom)
        assert result_tile.tx == tile.tx and \
            result_tile.ty == 0 and \
            result_zoom == zoom

    def test_tile_from_pixels(self):
        """Test making a tile from pixels"""
        gmp = GlobalMercatorProfile()
        zoom = 1
        point = MetersPoint(-20037508.3204, -20037508.3204)
        x, y = gmp.pixels_from_units(point, zoom)
