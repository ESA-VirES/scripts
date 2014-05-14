import sys
import os
from argparse import ArgumentParser
from itertools import izip
import csv
from itertools import tee, izip
from os.path import join, basename

import numpy as np
from spacepy import pycdf
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from lxml import etree

from osgeo import gdal

gdal.UseExceptions()
gdal.AllRegister()


report_template = """\
<?xml version='1.0' encoding='UTF-8'?>
<rep:browseReport xmlns:rep="http://ngeo.eo.esa.int/schema/browseReport" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://ngeo.eo.esa.int/schema/browseReport http://ngeo.eo.esa.int/schema/browseReport/browseReport.xsd" version="2.0">
  <rep:responsibleOrgName>EOX</rep:responsibleOrgName>
  <rep:dateTime>2014-01-24T14:53:33Z</rep:dateTime>
  <rep:browseType>SWARM-L1b</rep:browseType>
  <rep:browse>
    <rep:browseIdentifier/>
    <rep:fileName>%(browse_filename)s</rep:fileName>
    <rep:imageType>TIFF</rep:imageType>
    <rep:referenceSystemIdentifier>EPSG:4326</rep:referenceSystemIdentifier>
    <rep:verticalCurtainFootprint nodeNumber="%(node_number)s">
      <rep:colRowList>%(col_row_list)s</rep:colRowList>
      <rep:coordList>%(coord_list)s</rep:coordList>
      <rep:lookAngle>0.0</rep:lookAngle>
      <rep:verticalCurtainReferenceGrid>
        <rep:levelsNumber>104</rep:levelsNumber>
        <rep:heightLevelsList>239.232 479.056 718.893 958.726 1198.527 1438.38 1678.198 1918.025 2157.863 2397.685 2637.507 2877.35 3117.17 3357.011 3596.836 3836.667 4076.491 4316.328 4556.144 4795.979 5035.815 5275.637 5515.458 5755.296 5995.117 6234.956 6474.783 6714.608 6954.445 7194.282 7434.098 7673.934 7913.766 8153.581 8393.422 8633.252 8873.07 9112.904 9352.737 9592.547 9832.395 10072.215 10312.042 10551.875 10791.702 11031.519 11271.369 11511.196 11751.025 11990.863 12230.685 12470.507 12710.341 12950.161 13189.998 13429.829 13669.651 13909.474 14149.311 14389.13 14628.969 14868.805 15108.636 15348.458 15588.296 15828.117 16067.956 16307.781 16547.602 16787.434 17027.271 17267.081 17506.915 17746.753 17986.572 18226.408 18466.243 18706.068 18945.904 19185.737 19425.547 19665.395 19905.215 20145.042 20384.875 20624.702 20864.516 21104.358 21344.178 21584.012 21823.836 22063.667 22303.491 22543.328 22783.152 23022.989 23262.828 23502.651 23742.474 23982.311 24222.13 24461.967 24701.793 24941.618</rep:heightLevelsList>
      </rep:verticalCurtainReferenceGrid>
    </rep:verticalCurtainFootprint>
    <rep:startTime>%(start_time)s</rep:startTime>
    <rep:endTime>%(end_time)s</rep:endTime>
  </rep:browse>
</rep:browseReport>
"""

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


def savefig_pix(fig, fname, width, height, dpi=100):
    rdpi = 1.0/float(dpi)
    fig.set_size_inches(width*rdpi, height*rdpi)
    fig.savefig(fname, dpi=dpi)

def to_array(fig, width, height, dpi=100):
    rdpi = 1.0/float(dpi)
    fig.set_size_inches(width*rdpi, height*rdpi)
    fig.canvas.draw()
    data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    return data.reshape(fig.canvas.get_width_height()[::-1] + (3,))


def generate(input_filename, output_template, 
             start=0, stop=None, step=1000, x_size=1000, y_size=104):
    cdf = pycdf.CDF(input_filename)
    stop = stop or len(cdf["Timestamp"])

    r = range(start, stop, step)
    for i, (chunk_start, chunk_end) in enumerate(pairwise(r)):
        time = cdf["Timestamp"][chunk_start:chunk_end]
        lons = cdf["Longitude"][chunk_start:chunk_end]
        lats = cdf["Latitude"][chunk_start:chunk_end]
        rad = cdf["Radius"][chunk_start:chunk_end]
        f = cdf["F"][chunk_start:chunk_end]


        browse_filename = output_template % i
        report_filename = browse_filename.rpartition(".")[0] + ".xml"
        tmp_filename = browse_filename.rpartition(".")[0] + "_tmp_." + browse_filename.rpartition(".")[-1]
        
        fig = plt.figure()
        ax = fig.add_subplot(111, axisbg="black")

        x = np.arange(0, len(time))
        ax.fill_between(x, 0, f, facecolor="white")
        ax.plot(x, f, color="white")

        ax.set_xlim([0, len(time)])

        fig.subplots_adjust(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)
        savefig_pix(fig, tmp_filename, x_size, y_size)
        plt.close()

        # save TIFF via GDAL
        tmp_ds = gdal.Open(tmp_filename)
        tmp_band = tmp_ds.GetRasterBand(1)
        driver = gdal.GetDriverByName("GTiff")
        ds = driver.Create(browse_filename, x_size, y_size, 1, gdal.GDT_Byte)
        band = ds.GetRasterBand(1)
        band.WriteArray(tmp_band.ReadAsArray())
        del band, ds, tmp_band, tmp_ds

        os.remove(tmp_filename)

        # create a report

        out_step = step / 100

        # generate coord list
        
        col_row_list = " ".join(
            "%d 0" % i for i in xrange(0, step, out_step)
        )
        col_row_list += " %d 0" % (step - 1)

        coord_list = " ".join(
            "%f %f" % (lat, lon)
            for lat, lon in izip(lats[::out_step], lons[::out_step])
        )
        coord_list += " %f %f" % (lats[-1], lons[-1]) # finish with exact end


        node_number = len(coord_list.split(" ")) / 2


        report = report_template % {
            "node_number": node_number,
            "browse_filename": browse_filename,
            "col_row_list": col_row_list,
            "coord_list": coord_list,
            "start_time": time[0].isoformat("T"),
            "end_time": time[-1].isoformat("T"),
        }

        with open(report_filename, "w") as f:
            f.write(report)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--dir", "-d", default=".", help="Output directory")
    parser.add_argument("--offset", "-o", type=int, default=0, help="Start offset.")
    parser.add_argument("--stop", "-s", type=int, default=None, help="Stop offset")
    parser.add_argument("--chunksize", "-c", type=int, default=1000, help="Chunksize")
    parser.add_argument("--width", type=int, default=1000)
    parser.add_argument("--height", type=int, default=104)

    parser.add_argument("filename", nargs=1)
    parsed = parser.parse_args()
    generate(
        parsed.filename[0], join(parsed.dir, basename(parsed.filename[0] + "_%d.tif")), 
        parsed.offset, parsed.stop, parsed.chunksize,
        parsed.width, parsed.height
    )

