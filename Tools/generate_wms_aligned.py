#!/usr/bin/python

# Author: Steven D. Lander, RGi
# March 18, 2016

# This script assumes world-referenced tile coordinates

from collections import namedtuple

TileRange = namedtuple("TileRange", ["min", "max"])

def tile_bounds_geodetic(z, x, y):
    """Shamelessly taken from gdal2tiles.py"""
    tile_size = 256
    res_fact = 360.0 / tile_size
    res = res_fact / 2**z
    def calc(axis, max): return axis * tile_size * res - max
    return (calc(x, 180), calc(y, 90),
            calc((x + 1), 180), calc((y + 1), 90))

def iterate_tiles(z, range_x, range_y, **kwargs):
    try:
        base_url = kwargs["base_url"]
    except(AttributeError, KeyError):
        base_url = ""
    if type(range_x) is not TileRange or \
            type(range_y) is not TileRange:
        raise KeyError("z/x/y ranges must be TileRange objects")
    tmpl = "{}&BBOX={},{},{},{}"
    tile_urls = []
    for x in range(range_x.min, range_x.max+1):
        for y in range(range_y.min, range_y.max+1):
            # convert tile into bbox
            bbox = tile_bounds_geodetic(z, x, y)
            tile_urls.append(tmpl.format(base_url,
                bbox[0], bbox[1], bbox[2], bbox[3]))
    return tile_urls

# Print all these URLs out to the console, so they can
# be piped to output if desired
base_url = "http://localhost/GPEP/Hybrid-Performance-Test/service?"
base_url += "VERSION=1.3.0&REQUEST=GetMap&CRS=CRS:84&WIDTH=256&HEIGHT=256"
base_url += "&LAYERS=2,6,10,11,12&STYLES=,,,,&EXCEPTIONS=xml&FORMAT=image/jpeg"
base_url += "&BGCOLOR=0xFEFFFF&TRANSPARENT=TRUE"
tasks = [
        (6, TileRange(32, 55), TileRange(16,27)),
        (7, TileRange(72, 83), TileRange(40, 43)),
        (8, TileRange(160, 175), TileRange(64, 67)),
        (9, TileRange(336, 355), TileRange(128, 139)),
        (10, TileRange(676, 703), TileRange(256, 271)),
        (11, TileRange(1364, 1387), TileRange(524, 539)),
        (12, TileRange(2744, 2751), TileRange(1056, 1067)),
        (13, TileRange(5496, 5519), TileRange(2120, 2127)),
        (14, TileRange(11012, 11027), TileRange(4244, 4255)),
        (16, TileRange(44088, 44119), TileRange(17012, 17015)),
        ]

for task in tasks:
    level = iterate_tiles(task[0], task[1], task[2], base_url=base_url)
    for entry in level:
        print(entry)
