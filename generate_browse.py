import sys
from argparse import ArgumentParser
from itertools import izip
import csv

import numpy as np
from spacepy import pycdf
from matplotlib import pyplot as plt
from lxml import etree




def generate(input_filename, chunksize=1000, start=None, end=None, step=None,
             x_size=1000, y_size=100):
    cdf = pycdf.CDF(input_filename)

    time = cdf["Timestamp"][start:end]
    lon = cdf["Longitude"][start:end]
    lat = cdf["Latitude"][start:end]
    rad = cdf["Radius"][start:end]
    f = cdf["F"][start:end]


    fig = plt.figure()
    ax = fig.add_subplot(111, axisbg="black")
    
    x = np.arange(0, len(time))
    ax.plot(x, f)
    ax.set_xlim([0, len(time)])
    plt.axis("off")
    plt.savefig("out.png", bbox_inches='tight',) #facecolor="black")
    #plt.show()





if __name__ == "__main__":

    #parser = ArgumentParser()
    #parser.add_argument("filename", nargs=)
    generate(sys.argv[1], 10000, 0, 10000)
