#!/usr/bin/python

# Author: Steven D. Lander, RGi
# March 18, 2016

# This script assumes world-referenced tile coordinates

def gen_bbox(bbox, tiles):
    if tiles == 1:
        return ["&BBOX={},{},{},{}".format(bbox[0],
            bbox[1], bbox[2], bbox[3])]
    else:
        sub_count = tiles / 4
        min_x, min_y, max_x, max_y = bbox[0], bbox[1], bbox[2], bbox[3]
        middle_x = ((max_x - min_x) / 2) + min_x
        middle_y = ((max_y - min_y) / 2) + min_y
        bbox_ul = [min_x,    middle_y, middle_x, max_y]
        bbox_ur = [middle_x, middle_y, max_x, max_y]
        bbox_ll = [min_x,    min_y, middle_x, middle_y]
        bbox_lr = [middle_x, min_y, max_x, middle_y]
        return gen_bbox(bbox_ul, sub_count) + \
                gen_bbox(bbox_ur, sub_count) + \
                gen_bbox(bbox_ll, sub_count) + \
                gen_bbox(bbox_lr, sub_count)

base_url = "http://localhost/GPEP/Hybrid-Performance-Test/service?"
base_url += "VERSION=1.3.0&REQUEST=GetMap&CRS=CRS:84&WIDTH=256&HEIGHT=256"
base_url += "&LAYERS=2&STYLES=,,,,&EXCEPTIONS=xml&FORMAT=image/png"
base_url += "&BGCOLOR=0xFEFFFF&TRANSPARENT=TRUE"
z_min = 0
z_max = 7
tasks = [(2**x) * (2**x) for x in range(z_min, z_max+1)]
bbox = [50.92, 20.63, 78.12, 41.62]

bbox_list = []
for task in tasks:
    for entry in gen_bbox(bbox, task):
        print(base_url + entry)
