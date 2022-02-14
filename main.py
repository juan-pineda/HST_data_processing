# -- coding: utf-8 --
# Author: Juan Carlos Basto Pineda

import numpy as np
import astropy.convolution
from astropy.io import fits
import matplotlib.pyplot as plt

from functions import *

# 1-sigma PSF stimated from gaussian fitting to some stars in each band
# This were estimted for 60mas images and are given in pixel
new_sigma_hst =  {"105":1.47,
              "125":1.49,
              "140":1.55,
              "160":1.56,
              "435":0.78,
              "606":0.95,
              "775":0.87,
              "814":0.91,
              "850":0.85}


# Estimate quadratic difference to equate resolutions
def sigma_quad_diff(band):
    sigma = np.sqrt(new_sigma_hst['160']**2 - new_sigma_hst[band]**2)
    return sigma


# Read an image and creates the low resolution version, 
# sky-subtracted and masked. 
# There is the option to convolve with the gaussian kernel squared
def signal_low_resolution(gal,band,case='60', kernel_square=False):
    hdul = new_read_hst(gal,band,case)
    sky = read_sky(gal,band)
    new_mask = enlarge_mask(gal)
    sigma = sigma_quad_diff(band)
    data = hdul[0].data - sky
    data = data * new_mask
    data = direct_convolution(data,sigma,kernel_square)
    data = data * new_mask
    return data


# Read a variance map and creates the low resolution version, masked
def variance_low_resolution(gal,band,case='60'):
    weight = read_hst_weight_maps(gal,band,case)
    sigma = sigma_quad_diff(band)
    variance = 1./weight[0].data
    variance = direct_convolution(variance,sigma,kernel_square=True)
    new_mask = enlarge_mask(gal)
    variance = variance*new_mask
    return variance


def direct_convolution(data,sigma,kernel_square=False):
    FWHM = sigma * (2 * np.sqrt(2 * np.log(2)))
    FWHM = np.int(FWHM)+1
    data2 = np.zeros((data.shape[0] + 2 * FWHM, data.shape[1] + 2 * FWHM))
    data2[FWHM:data.shape[0] + FWHM, FWHM:data.shape[1] + FWHM] = data
    y, x = np.indices((data2.shape))
    psf = (1. / (2 * np.pi * sigma ** 2)) * np.exp(-((x - data2.shape[1] / 2) ** 2 + (y - data2.shape[0] / 2) ** 2) / (2.0 * sigma ** 2))
    psf /= psf.sum()  # normalisation PSF
    if kernel_square:
        psf = psf**2
    psf_shift = fftshift(psf)
    data_conv = irfft2(rfft2(data2) * rfft2(psf_shift))
    data_conv = data_conv[FWHM:data.shape[0] + FWHM, FWHM:data.shape[1] + FWHM]
    return data_conv






