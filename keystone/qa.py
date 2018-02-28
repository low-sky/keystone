from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os

def plotCoords(filename, plotdir=None):
    if not plotdir:
        plotdir = os.getcwd()

    hdulist = fits.open(filename)
    s = hdulist[1]
    lon = s.data['CRVAL2']
    lat = s.data['CRVAL3']
    fig =  plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.plot(lon, lat,'ro')
    ax.set_xlabel('Tracking System Longitude (deg)')
    ax.set_ylabel('Tracking System Latitude (deg)')
    ax.set_title(filename)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(plotdir + '/' + filename.replace('.fits','_encoder.png'))
    plt.close()
    plt.clf()
    
