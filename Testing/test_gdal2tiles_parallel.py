#!/usr/bin/python

from sys import path
from os.path import abspath
path.append(abspath("Tiling"))

from gdal2tiles_parallel import Tile
from gdal2tiles_parallel import MetersPoint
from gdal2tiles_parallel import PixelsPoint
from gdal2tiles_parallel import LonLatPoint
from gdal2tiles_parallel import MetersToTile
from gdal2tiles_parallel import ITileProfile
from gdal2tiles_parallel import GlobalMercatorProfile

class TestITileProfile:

    def test_lower_left_tile(self):
        itp = ITileProfile()
        tile = Tile(0,0)
        zoom = 1
        result = itp.lower_left_tile(tile, zoom)
        assert result.tx == 0 and result.ty == 0

    def test_upper_left_tile_1(self):
        itp = ITileProfile()
        tile = Tile(0, 0)
        zoom = 1
        result = itp.upper_left_tile(tile, zoom)
        assert result.tx == 0 and result.ty == 1

    def test_upper_left_tile_2(self):
        itp = ITileProfile()
        tile = Tile(0, 31)
        zoom = 13
        result = itp.upper_left_tile(tile, zoom)
        assert result.tx == 0 and result.ty == 8160

    def test_quad_tree(self):
        itp = ITileProfile()
        tile = Tile(0, 31)
        zoom = 13
        result = itp.quad_tree(tile, zoom)
        assert result == '2222222200000'

class TestGlobalMercatorProfile:

    def test_init(self):
        gmp = GlobalMercatorProfile()
        assert gmp.tile_size == 256 and \
            gmp.initial_resolution == 156543.03392804097 and \
            gmp.origin_shift == 20037508.342789244

    def test_resolution(self):
        gmp = GlobalMercatorProfile()
        result = gmp.resolution(13)
        assert result == 19.109257071294063

    def test_zoom_for_pixel_size(self):
        gmp = GlobalMercatorProfile()
        pixel_size = 10
        result = gmp.zoom_for_pixel_size(pixel_size)
        assert result == 13

    def test_pixels_to_tile_1(self):
        gmp = GlobalMercatorProfile()
        point = PixelsPoint(256, 256)
        result = gmp.to_tile(point)
        assert result.tx == 0 and result.ty == 0

    def test_pixels_to_tile_2(self):
        gmp = GlobalMercatorProfile()
        point = PixelsPoint(256, 257)
        result = gmp.to_tile(point)
        assert result.tx == 0 and result.ty == 1

    def test_meters_to_tile_1(self):
        gmp = GlobalMercatorProfile()
        point = MetersPoint(-20037508.34, 20037508.34)
        result = gmp.to_tile(point, zoom=0)
        assert result.tx == 0 and result.ty == 0

    def test_meters_to_tile_2(self):
        gmp = GlobalMercatorProfile()
        point = MetersPoint(7682838.58, 4106808.65)
        result = gmp.to_tile(point, zoom=14)
        print(result.tx, result.ty)
        assert result.tx == 11332 and result.ty == 9870

    def test_tile_to_tile(self):
        gmp = GlobalMercatorProfile()
        tile = Tile(34, 100)
        result = gmp.to_tile(tile)
        assert result.tx == tile.tx and \
            result.ty == tile.ty

    def test_lon_lat_to_tile(self):
        gmp = GlobalMercatorProfile()
        point = LonLatPoint(-34.123, 78.1234)
        try:
            gmp.to_tile(point)
            assert False
        except NotImplementedError as e:
            assert e is not None and type(e) == NotImplementedError
