# -- coding: utf-8 --
# Author: Juan Carlos Basto Pineda

""" Config params and main functions """

import numpy as np
import astropy.convolution
from astropy.io import fits
import matplotlib.pyplot as plt

from functions import *

# 1-sigma PSF stimated from gaussian fitting to some stars in each band
# this values were estimated from the 60mas images and are given in pixels
sigma_hst =  {"105":1.47,
              "125":1.49,
              "140":1.55,
              "160":1.56,
              "435":0.78,
              "606":0.95,
              "775":0.87,
              "814":0.91,
              "850":0.85}


# quadratic difference between kernels to equate resolutions
# 160 is the poorest-resolution band
def sigma_quad_diff(band):
    sigma = np.sqrt(sigma_hst['160']**2 - sigma_hst[band]**2)
    return sigma


# create the low-resolution / masked / sky-subtracted version of an image
def signal_low_resolution(gal, band):
    hdul = new_read_hst(gal, band, case='60')
    sky = read_sky(gal,band)
    flux = hdul[0].data - sky
    # we expand the mask a little bit to not miss pixels in the border
    new_mask = enlarge_mask(gal)
    # and multiply to avoid sky contamination in the convolution
    flux = flux * new_mask 

    sigma = sigma_quad_diff(band)
    flux = direct_convolution(flux, sigma, kernel_square=False)
    flux = flux * new_mask
    return flux


# read a variance map and creates the low-resolution / masked version
def variance_low_resolution(gal, band):
    weight = read_hst_weight_maps(gal, band, case='60')
    variance = 1./weight[0].data

    sigma = sigma_quad_diff(band)
    variance = direct_convolution(variance, sigma, kernel_square=True)
    new_mask = enlarge_mask(gal)
    variance = variance*new_mask
    return variance


# sigma must be in pixels
# use the kernel squared for the variance
def direct_convolution(data, sigma, kernel_square=False):
    FWHM = sigma * (2 * np.sqrt(2 * np.log(2)))
    FWHM = np.int(FWHM)+1
    data2 = np.zeros((data.shape[0] + 2 * FWHM, data.shape[1] + 2 * FWHM))
    data2[FWHM:data.shape[0] + FWHM, FWHM:data.shape[1] + FWHM] = data
    y, x = np.indices((data2.shape))

    # gaussian kernel
    exponent = -((x - data2.shape[1] / 2) ** 2 + (y - data2.shape[0] / 2) ** 2)
    exponent = exponent / (2.0 * sigma ** 2)
    psf = (1. / (2 * np.pi * sigma ** 2)) * np.exp(exponent)
    # normalisation of the PSF
    psf /= psf.sum()

    if kernel_square:
        psf = psf**2
        
    psf_shift = fftshift(psf)
    data_conv = irfft2(rfft2(data2) * rfft2(psf_shift))
    data_conv = data_conv[FWHM:data.shape[0] + FWHM, FWHM:data.shape[1] + FWHM]
    return data_conv






