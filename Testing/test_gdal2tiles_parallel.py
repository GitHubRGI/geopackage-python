#!/usr/bin/python

from sys import path
from os.path import abspath
path.append(abspath("Tiling"))

from gdal2tiles_parallel import Tile
from gdal2tiles_parallel import MetersPoint
from gdal2tiles_parallel import PixelsPoint
from gdal2tiles_parallel import LonLatPoint
from gdal2tiles_parallel import ITileProfile
from gdal2tiles_parallel import GlobalMercatorProfile

class TestITileProfile:

    def test_lower_left_tile_1(self):
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

    def test_quad_tree_1(self):
        itp = ITileProfile()
        tile = Tile(0, 31)
        zoom = 13
        result = itp.quad_tree(tile, zoom)
        assert result == '2222222200000'

class TestGlobalMercatorProfile:

    def test_init_1(self):
        gmp = GlobalMercatorProfile()
        assert gmp.tile_size == 256 and \
            gmp.initial_resolution == 156543.03392804097 and \
            gmp.origin_shift == 20037508.342789244

    def test_init_2(self):
        gmp = GlobalMercatorProfile(tile_size=512)
        assert gmp.tile_size == 512 and \
            gmp.initial_resolution == 78271.51696402048

    def test_resolution_1(self):
        gmp = GlobalMercatorProfile()
        result = gmp.resolution(13)
        assert result == 19.109257071294063

    def test_zoom_for_pixel_size_1(self):
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
        assert result.tx == 11332 and result.ty == 9870

    def test_meters_to_tile_3(self):
        gmp = GlobalMercatorProfile()
        point = MetersPoint(7682838.58, 4106808.65)
        try:
            gmp.to_tile(point)
            assert False
        except KeyError as e:
            assert e is not None and type(e) == KeyError

    def test_meters_to_tile_4(self):
        gmp = GlobalMercatorProfile()
        point = MetersPoint(7682838.58, 4106808.65)
        try:
            gmp.to_tile(point, foo=12)
            assert False
        except KeyError as e:
            assert e is not None and type(e) == KeyError

    def test_tile_to_tile_1(self):
        gmp = GlobalMercatorProfile()
        tile = Tile(34, 100)
        result = gmp.to_tile(tile)
        assert result.tx == tile.tx and \
            result.ty == tile.ty

    def test_lon_lat_to_tile_1(self):
        gmp = GlobalMercatorProfile()
        point = LonLatPoint(-34.123, 78.1234)
        try:
            gmp.to_tile(point)
            assert False
        except NotImplementedError as e:
            assert e is not None and type(e) == NotImplementedError

    def test_lon_lat_to_map_coords_1(self):
        gmp = GlobalMercatorProfile()
        point = LonLatPoint(180, 85)
        result = gmp.to_map_coordinates(point)
        assert result.mx == 20037508.342789244 and \
            result.my == 19971868.88040853

    def test_lon_lat_to_map_coords_2(self):
        gmp = GlobalMercatorProfile()
        point = LonLatPoint(-148.99, 61.60)
        result = gmp.to_map_coordinates(point)
        assert result.mx == -16585490.933289833 and \
            result.my == 8764912.67945232

    def test_pixels_to_map_coords_1(self):
        gmp = GlobalMercatorProfile()
        point = PixelsPoint(256, 256)
        result = gmp.to_map_coordinates(point, zoom=0)
        assert result.mx == 20037508.342789244 and \
            result.my == 20037508.342789244

    def test_pixels_to_map_coords_2(self):
        gmp = GlobalMercatorProfile()
        point = PixelsPoint(0, 0)
        result = gmp.to_map_coordinates(point, zoom=0)
        assert result.mx == -20037508.342789244 and \
            result.my == -20037508.342789244

    def test_pixels_to_map_coords_3(self):
        gmp = GlobalMercatorProfile()
        point = PixelsPoint(0, 0)
        try:
            gmp.to_map_coordinates(point)
            assert False
        except KeyError as e:
            assert e is not None and type(e) == KeyError

    def test_pixels_to_map_coords_4(self):
        gmp = GlobalMercatorProfile()
        point = PixelsPoint(0, 0)
        try:
            gmp.to_map_coordinates(point, foo=12)
            assert False
        except KeyError as e:
            assert e is not None and type(e) == KeyError

    def test_meters_to_map_coords_1(self):
        gmp = GlobalMercatorProfile()
        point = MetersPoint(41, 182)
        result = gmp.to_map_coordinates(point)
        assert point.mx == result.mx and \
            point.my == result.my

    def test_tile_to_map_coordinates_1(self):
        gmp = GlobalMercatorProfile()
        tile = Tile(1, 2)
        try:
            gmp.to_map_coordinates(tile)
            assert False
        except NotImplementedError as e:
            assert e is not None and type(e) == NotImplementedError

    def test_meters_to_lon_lat_1(self):
        gmp = GlobalMercatorProfile()
        point = MetersPoint(20037508.342789244, 20037508.342789244)
        result = gmp.to_lon_lat(point)
        assert result.lon == 180.0 and \
            result.lat == 85.0511287798066

    def test_lon_lat_to_lon_lat_1(self):
        gmp = GlobalMercatorProfile()
        point = LonLatPoint(123.4567890, 98.7654321)
        result = gmp.to_lon_lat(point)
        assert point.lon == result.lon and \
            point.lat == result.lat

    def test_pixels_to_lon_lat_1(self):
        gmp = GlobalMercatorProfile()
        point = PixelsPoint(1, 2)
        try:
            gmp.to_lon_lat(point)
            assert False
        except NotImplementedError as e:
            assert e is not None and type(e) == NotImplementedError

    def test_meters_to_pixels_1(self):
        gmp = GlobalMercatorProfile()
        point = MetersPoint(20037508.342789244, 20037508.342789244)
        result = gmp.to_pixels(point, zoom=0)
        assert result.x == 256 and result.y == 256

    def test_meters_to_pixels_2(self):
        gmp = GlobalMercatorProfile()
        point = MetersPoint(-20037508.342789244, -20037508.342789244)
        result = gmp.to_pixels(point, zoom=0)
        assert result.x == 0 and result.y == 0

    def test_meters_to_pixels_3(self):
        gmp = GlobalMercatorProfile()
        point = MetersPoint(-20037508.342789244, -20037508.342789244)
        try:
            gmp.to_pixels(point)
            assert False
        except KeyError as e:
            assert e is not None and type(e) == KeyError

    def test_meters_to_pixels_4(self):
        gmp = GlobalMercatorProfile()
        point = MetersPoint(-20037508.342789244, -20037508.342789244)
        try:
            gmp.to_pixels(point, foo=12)
            assert False
        except KeyError as e:
            assert e is not None and type(e) == KeyError

    def test_pixels_to_pixels_1(self):
        gmp = GlobalMercatorProfile()
        point = PixelsPoint(1, 2)
        result = gmp.to_pixels(point, zoom=12)
        assert point.x == result.x and point.y == result.y

    def test_lon_lat_to_pixels_1(self):
        gmp = GlobalMercatorProfile()
        point = LonLatPoint(123.4567890, 98.7654321)
        try:
            gmp.to_pixels(point, zoom=18)
            assert False
        except NotImplementedError as e:
            assert e is not None and type(e) == NotImplementedError

    def test_pixels_to_raster_1(self):
        gmp = GlobalMercatorProfile()
        point = PixelsPoint(256, 256)
        [x, y] = gmp.to_raster(point, zoom=0)
        assert x == 256 and y == 0

    def test_pixels_to_raster_2(self):
        gmp = GlobalMercatorProfile()
        point = PixelsPoint(512, 0)
        [x, y] = gmp.to_raster(point, zoom=1)
        assert x == 512 and y == 512

    def test_pixels_to_raster_3(self):
        gmp = GlobalMercatorProfile()
        point = PixelsPoint(512, 0)
        try:
            gmp.to_raster(point, foo=1)
            assert False
        except KeyError as e:
            assert e is not None and type(e) == KeyError

    def test_lon_lat_to_raster_1(self):
        gmp = GlobalMercatorProfile()
        point = LonLatPoint(180.0, 90.0)
        try:
            gmp.to_raster(point, zoom=0)
            assert False
        except NotImplementedError as e:
            assert e is not None and type(e) == NotImplementedError
