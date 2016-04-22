#!/usr/bin/python

# Author: Steven D. Lander, RGi
# March 17, 2016

from collections import namedtuple

TileRange = namedtuple("TileRange", ["min", "max"])

def iterate_tiles(z, range_x, range_y, **kwargs):
    try:
        base_url = kwargs["base_url"]
    except(AttributeError, KeyError):
        base_url = ""
    try:
        file_ext = kwargs["file_ext"]
        if "." not in file_ext:
            raise ValueError("file_ext needs a period in it.")
    except(AttributeError, KeyError):
        file_ext = ""
    except(ValueError):
        raise
    if type(range_x) is not TileRange or \
            type(range_y) is not TileRange:
        raise KeyError("z/x/y ranges must be TileRange objects")
    tmpl = "{}{}/{}/{}{}"
    tile_urls = []
    for x in range(range_x.min, range_x.max+1):
        for y in range(range_y.min, range_y.max+1):
            tile_urls.append(tmpl.format(base_url, z, x, y, file_ext))
    return tile_urls

# Print all these URLs out to the console, so they can
# be piped to output if desired
base_url = "http://example/url/"
file_ext = ".png"
tasks = [
        (1, TileRange(0, 3), TileRange(0, 3)), # All zoom 1 tiles (EPSG:3857)
        # (7, TileRange(0, 2**7), TileRange(0, 2**(7-1))), # All zoom 7 tiles (EPSG:4326)
        # (16, TileRange(0, 2**16), TileRange(0, 2**16)), # All zoom 16 tiles (EPSG:3857)
        ]

for task in tasks:
    level = iterate_tiles(task[0], task[1], task[2], base_url=base_url,
            file_ext=file_ext)
    for entry in level:
        print(entry)
