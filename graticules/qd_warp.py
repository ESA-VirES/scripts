#!/usr/bin/env python
"""Computes magnetic field graticules."""

# Authors:  Martin Paces <martin.paces@eox.at>
#           Joachim Ungar <joachim.ungar@eox.at>
#
# ------------------------------------------------------------------------------
# Copyright (C) 2016 EOX IT Services GmbH
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies of this Software or works derived from this Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ------------------------------------------------------------------------------

import sys
import os
import fiona
import argparse
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, mapping, shape
from affine import Affine
import rasterstats
from eoxmagmod import convert, GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL
from eoxmagmod.qd import eval_qdlatlon


def main(args):
    """Parse arguments and run graticule creation."""
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=str, help="output vector data")
    parser.add_argument("stepsize", type=int, help="degree step size")
    parser.add_argument(
        "--size_x", type=int, help="model raster x size", default=8192)
    parser.add_argument(
        "--size_y", type=int, help="model raster y size", default=4096)
    parser.add_argument(
        "--driver", type=str, help="output driver", default="ESRI Shapefile")
    parsed = parser.parse_args(args)

    output_file = parsed.output
    stepsize = parsed.stepsize
    size_x = parsed.size_x
    size_y = parsed.size_y

    fieldname = "dd"
    elevation = 0
    bounds = (-180., -90., 180., 90.)
    left, bottom, right, top = bounds
    pixel_x_size = (right - left) / float(size_x)
    pixel_y_size = (top - bottom) / float(size_y)
    affine = Affine.translation(left, top) * Affine.scale(
        pixel_x_size, -pixel_y_size)

    if os.path.isfile(output_file):
        os.remove(output_file)

    lons, lats = np.meshgrid(
        np.linspace(left, right, size_x, endpoint=True),
        np.linspace(top, bottom, size_y, endpoint=True)
    )

    # Geodetic coordinates with elevation above the WGS84 ellipsoid.
    coord_gdt = np.empty((size_y, size_x, 3))
    coord_gdt[:, :, 0] = lats
    coord_gdt[:, :, 1] = lons
    coord_gdt[:, :, 2] = elevation
    decimal_year = np.empty((size_y, size_x))
    decimal_year.fill(2016.1)
    coord_gct = convert(coord_gdt, GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL)
    qd_lat, qd_lon = eval_qdlatlon(
        coord_gct[..., 0].ravel(), coord_gct[..., 1].ravel(),
        coord_gct[..., 2].ravel(), decimal_year.ravel())

    # Latitudes
    lat = np.reshape(qd_lat, (size_y, size_x))
    lat_lines = extract_contours(lat, bounds, stepsize, fieldname, 0)

    # Longitudes
    lon = np.reshape(qd_lon, (size_y, size_x))
    lon_lines = extract_contours(
        np.absolute(lon), bounds, stepsize, fieldname, 0)
    meridians = extract_contours(lon, bounds, stepsize, fieldname, 0)
    # for row in range(len(lon)):
    #     diff = lon[row][:-1]-lon[row][1:]
    #     step_idxes = np.argwhere(diff > 180.)
    #     if step_idxes:
    #         lon[row][np.argwhere(diff > 180.)[0]+1:] += 360.
    # lon_lines = extract_contours(lon, bounds, stepsize, fieldname)

    out_schema = dict(
        geometry="LineString",
        properties=dict(
            degrees="int", direction="str", display="str", dd="float",
            scalerank="int"))
    with fiona.open(
        output_file, "w", schema=out_schema, crs={'init': 'epsg:4326'},
        driver=parsed.driver
    ) as dst:
        for feature in lat_lines:
            latitude = feature["properties"]["dd"]
            lat_int = int(round(latitude))
            direction = "N" if lat_int > 0 else "S"
            display_lat = lat_int if lat_int > 0 else -lat_int
            display = "%s %s" % (display_lat, direction)
            if lat_int == 0:
                direction = None
                display = "0"
            feature["properties"].update(
                degrees=lat_int, direction=direction, display=display,
                scalerank=None)
            dst.write(feature)
        for feature in _extract_longitudes(lon_lines, lon, affine):
            dst.write(feature)
            longitude = feature["properties"]["dd"]
            lon_int = int(round(longitude))
            direction = "W" if lon_int > 0 else "E"
            display_lon = lon_int if lon_int > 0 else -lon_int
            display = "%s %s" % (display_lon, direction)
            if lon_int == 0:
                direction = None
                display = "0"
            feature["properties"].update(
                degrees=lon_int, direction=direction, display=display,
                scalerank=None)
            # if longitude > 180.:
            #     feature["properties"].update(dd=longitude-360.)
            # dst.write(feature)

        for feature in meridians:
            longitude = feature["properties"]["dd"]
            if longitude in [0, 180]:
                feature["properties"].update(
                    degrees=0, direction=None, display="0", scalerank=None)
                dst.write(feature)


def _extract_longitudes(lon_lines, lon_array, affine):
    for line in lon_lines:
        out_line_coords = []
        longitude = line["properties"]["dd"]
        for x, y, in shape(line["geometry"]).coords:
            point = Point(x, y)
            rastervalue = rasterstats.gen_point_query(
                point, lon_array, affine=affine).next()
            if rastervalue is None:
                continue

            if not out_line_coords:
                out_line_coords.append(point)
                if rastervalue >= 0.:
                    current_is_positive = True
                else:
                    current_is_positive = False

            if rastervalue >= 0.:
                if current_is_positive:
                    out_line_coords.append(point)
                else:
                    out_line_coords.append(point)
                    out_line = line
                    out_line.update(
                        geometry=mapping(LineString(out_line_coords)))
                    # longitude = -longitude
                    lon_int = int(round(-longitude))
                    direction = "W" if lon_int > 0 else "E"
                    display_lon = lon_int if lon_int > 0 else -lon_int
                    display = "%s %s" % (display_lon, direction)
                    if lon_int == 0:
                        direction = None
                        display = "0"
                    out_line["properties"].update(
                        dd=-longitude, degrees=lon_int, direction=direction,
                        display=display, scalerank=None)
                    out_line_coords = []
                    yield out_line
            elif rastervalue < 0.:
                if current_is_positive:
                    out_line_coords.append(point)
                    out_line = line
                    out_line.update(
                        geometry=mapping(LineString(out_line_coords)))
                    lon_int = int(round(longitude))
                    direction = "W" if lon_int > 0 else "E"
                    display_lon = lon_int if lon_int > 0 else -lon_int
                    display = "%s %s" % (display_lon, direction)
                    if lon_int == 0:
                        direction = None
                        display = "0"
                    out_line["properties"].update(
                        dd=longitude, degrees=lon_int, direction=direction,
                        display=display, scalerank=None)
                    out_line_coords = []
                    yield out_line
                else:
                    out_line_coords.append(point)

        if out_line_coords:
            if not current_is_positive:
                longitude = -longitude
            out_line = line
            out_line.update(
                geometry=mapping(LineString(out_line_coords)))
            lon_int = int(round(longitude))
            direction = "W" if lon_int > 0 else "E"
            display_lon = lon_int if lon_int > 0 else -lon_int
            display = "%s %s" % (display_lon, direction)
            if lon_int == 0:
                direction = None
                display = "0"
            out_line["properties"].update(
                dd=longitude, degrees=lon_int, direction=direction,
                display=display, scalerank=None)
            out_line_coords = []
            yield out_line


def extract_contours(array, bounds, interval=100, field='elev', round_=0):
    """
    Extract contour lines from an array in a given interval.

    Returns GeoJSON-like dictionary.
    """
    levels = _get_contour_values(
        array.min(), array.max(), interval=interval)
    if not levels:
        return []
    try:
        contours = plt.contour(array, levels)
    except:
        raise
    index = 0
    out_contours = []
    left, bottom, right, top = bounds
    pixel_x_size = (right - left) / (array.shape[1] - 1)
    pixel_y_size = (top - bottom) / (array.shape[0] - 1)
    for level in range(len(contours.collections)):
        elevation = levels[index]
        index += 1
        paths = contours.collections[level].get_paths()
        for path in paths:
            out_coords = [
                (
                    left+(i[1]*pixel_x_size),
                    top-(i[0]*pixel_y_size),
                )
                for i in zip(path.vertices[:, 1], path.vertices[:, 0])
            ]
            if len(out_coords) >= 2:
                line = LineString(out_coords)
                out_contours.append({
                        'properties': {field: elevation},
                        'geometry': mapping(line)
                })

    return out_contours


def _stretch(array, target_min=0., target_max=180.):
    """Stretch array values to new minimum and maximum."""
    assert target_min < target_max
    new_range = array.max() - array.min()
    old_range = target_max - target_min
    array = np.clip(array, array.min(), array.max())[:]
    return (
        (
            (array.astype(float) - array.min()) * (  # normalize to min
                float(old_range) / float(new_range)
            )
        )
    ).astype(array.dtype)


def _get_contour_values(min_val, max_val, base=0, interval=100):
    """Return a list of values between min and max within an interval."""
    i = base
    out = []

    if min_val < base:
        while i >= min_val:
            i -= interval

    while i <= max_val:
        if i >= min_val:
            out.append(float(i))
        i += interval

    return out


if __name__ == "__main__":
    main(sys.argv[1:])
