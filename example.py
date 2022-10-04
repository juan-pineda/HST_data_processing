# -- coding: utf-8 --
# Author: Juan Carlos Basto Pineda

""" Create a set of homogeneous images and weight maps for a set of galaxies """

from main import *

# directory to store the data products generated here
outdir = './results_example'

# process data for several galaxies
galaxies = ["3","15","37","912","919","937","943","982","1002","m3","m943"]

# Galaxies 3 and 943 are mergers. To retain only the main object we pass the
# number of the galaxy (e.g. 3), but we use an extra 'm' in front (e.g. m3)
# if we want to create maps of the merging objects as a single thing.

for gal in galaxies:
    for band in ["105","125","140","435","606","775","814","850"]:

        # we need to read this just to propagate the main info in the header
        weight = read_hst_weight_maps(gal,band)

        # low resolution version of the 60mas image
        newdata = signal_low_resolution(gal, band)
        filename = gal+'_'+band+'_lowres.fits'
        fits.writeto(outdir+'/'+filename, newdata, weight[0].header, overwrite=True)

        # low resolution version of the 60mas variance map
        newdata = variance_low_resolution(gal, band)
        filename = gal+'_'+band+'_variance_lowres.fits'
        fits.writeto(outdir+'/'+filename, newdata, weight[0].header, overwrite=True)

# WARNING:
# This routine was developed for 60mas only. Mixing 30mas data would require
# several fixes: diferent PSF/SKY, variance correction, scale pixels size, etc.



