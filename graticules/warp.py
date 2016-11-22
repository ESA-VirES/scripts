#!/usr/bin/env python
"""Converts input vector data into magnetically warped output vector data."""

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

import os
import sys
import argparse
import fiona
from shapely.geometry import (
    shape, mapping, box, LineString, MultiLineString, GeometryCollection)
from eoxmagmod import convert, GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL
from eoxmagmod.qd import eval_qdlatlon


def main(args):
    """Parse arguments and run warp."""
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str, help="input vector data")
    parser.add_argument("output", type=str, help="output vector data")
    parsed = parser.parse_args(args)

    if os.path.isfile(parsed.output):
        os.remove(parsed.output)

    with fiona.open(parsed.input, "r") as src:
        assert src.crs == {'init': u'epsg:4326'}
        with fiona.open(
            parsed.output, "w", schema=src.schema.copy(), driver=src.driver,
            crs=src.crs
        ) as dst:
            for feature in src:
                feature.update(
                    geometry=mapping(
                        warp_geometry(shape(feature["geometry"]))))
                dst.write(feature)


def warp_geometry(geom):
    """
    Warp geometry while preserving geometry type.

    Currently just works on LineString or MultiLineString geometries.
    """
    if geom.type == "LineString":
        return _warp_geometry(geom)
    elif geom.type == "MultiLineString":
        return MultiLineString([_warp_geometry(geom) for line in geom])
    else:
        raise IOError("invalid input geometry type: %s" % geom.type)


def _warp_geometry(geom):
    """Return warped geometry while correctly dealing with Antimeridian."""
    out_geom = LineString([_magnetic_warp(x, y) for x, y in geom.coords])
    # Geometries crossing the Antimeridian will get clipped and returned as
    # MultiLineStrings.
    if not WGS84_BOUNDS.contains(out_geom):
        out_geom = MultiLineString([
            line
            for line in [
                WGS84_BOUNDS.intersection(out_geom),
                _shift_geom(
                    WGS84_LEFT.intersection(out_geom), xoff=360, yoff=0),
                _shift_geom(
                    WGS84_RIGHT.intersection(out_geom), xoff=-360, yoff=0)]
            if line.type == "LineString"  # kick out empty GeometryCollection
        ])
    return out_geom


WGS84_BOUNDS = box(-180, -90, 180, 90)
WGS84_LEFT = box(-540, -90, -180, 90)
WGS84_RIGHT = box(180, -90, 540, 90)


def _magnetic_warp(x, y, decimal_year=2016.0):
    """Apply magnetic warp magic."""
    gc_lat, gc_lon, gc_rad = convert(
        [y, x, 0], GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL)
    qd_lat, qd_lon = eval_qdlatlon(gc_lat, gc_lon, gc_rad, decimal_year)
    return qd_lon, qd_lat


def _shift_geom(geom, xoff=0, yoff=0):
    """Shift geometry along x and/or y axis."""
    if geom.type == "LineString":
        return LineString([(x+xoff, y+yoff) for x, y in geom.coords])
    else:
        return GeometryCollection()


if __name__ == "__main__":
    main(sys.argv[1:])
