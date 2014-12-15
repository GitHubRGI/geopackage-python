geopackage-python : Python-based tools for creating OGC GeoPackages.
=================

[GeoPackage Specification from the Open Geospatial
Consortium](http://opengeospatial.org/standards/geopackage)


## gdal2tiles_parallel.py

[See PDF documentation
here.](https://github.com/GitHubRGI/geopackage-python/raw/master/Documentation/release/Instructions_For_gdal2tiles.pdf)


## tiles2gpkg_parallel.py

A script that will accept a folder full of tiles arranged in TMS or WMTS format
(z/x/y) and output an OGC-compliant GeoPackage.

### Requirements

- [Python 2.7](https://www.python.org/downloads)
- Python Imaging Library (PIL or Pillow) - `pip -U Pillow`

### Usage
[See PDF documentation
here.](https://github.com/GitHubRGI/geopackage-python/raw/master/Documentation/release/Instructions_For_tiles2gpkg.pdf)

### Testing

To run the test suite for tiles2gpkg_parallel.py, you will need pytest for
python 2.7. For detailed information on how to install or set up pytest, see
the [official documentation.](http://pytest.org/latest/getting-started.html)

Using PIP from the command line (Windows or Linux):

    pip install -U pytest

To run the test suite, copy tiles2gpkg_parallel.py into the Testing folder and
then run:

    py.test test_tiles2gpkg.py


## Known Issues

Refer to the issue tracker for updated information on bugs, planned features,
and other news.
