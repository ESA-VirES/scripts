import sys
from argparse import ArgumentParser
from itertools import izip
import csv

import numpy as np
from spacepy import pycdf
from matplotlib import pyplot as plt
from lxml import etree


def convert(input_filename, output_filename):
    cdf = pycdf.CDF(input_filename)
    values = izip(
        cdf["Timestamp"][...], cdf["Longitude"][...], cdf["Latitude"][...], cdf["Radius"][...], cdf["F"][...]
    )

    with open(output_filename, "w+") as f:
        writer = csv.writer(f)
        
        writer.writerow(("timestamp", "lon", "lat", "radius", "f"))
        for timestamp, lon, lat, radius, f in values:
            writer.writerow((timestamp.isoformat("T"), lon, lat, radius, f))


if __name__ == "__main__":
    convert(sys.argv[1], sys.argv[2])
