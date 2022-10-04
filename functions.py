# -- coding: utf-8 --
# Author: Juan Carlos Basto Pineda

""" Set paths / Deal with segmentation maps and create interpolated masks """ 

import os
import sys
import cv2
import pickle
import numpy as np
from astropy import wcs
import scipy.interpolate
from astropy.io import fits
from shutil import copyfile
from scipy import interpolate
import matplotlib.pyplot as plt
from numpy.fft import fftshift, rfft2, irfft2
from mpl_toolkits.axes_grid1 import make_axes_locatable


  ####--- Manual Configurations --------------####
  ####--- Set the paths to the source data ---####

new_folder_hst = "./INPUT_DATA/new_stamps/" # HST data
new_folder_sky = "./INPUT_DATA/NEW_SKY/" # measurements of the sky
segmendir = "./INPUT_DATA/segment_maps/" # segmentation masks

  ###--- End of manual configurations ------####
  ###---------------------------------------####


# The most important function here is ***enlarge_mask(gal)***
# It is called from "main.py" and makes use of all the other functions here
# Yet the code is presented in a sequence intended to ease its comprehension


# read HST image
def new_read_hst(gal, band, case='60'):
    """
    'case' stands for pixel size in mas, must be: '30' or '60'
    'm3' is the same galaxy 3 (see readme)
    """
    if gal == "m3":
        gal = "3"
    elif gal == "m943":
        gal = "943"
    try:
         hdul = fits.open(new_folder_hst+gal+'_'+case+'mas_f'+band+"w_sci.fits")
    except:
        hdul = fits.open(new_folder_hst+gal+'_'+case+'mas_f'+band+"lp_sci.fits")
    return hdul


# read weight map
def read_hst_weight_maps(gal, band, case='60'):
    """
    'case' stands for pixel size in mas, must be: '30' or '60'
    'm3' is the same galaxy 3 (see readme)
    """
    if gal == "m3":
        gal = "3"
    elif gal == "m943":
        gal = "943"
    try:
         hdul = fits.open(new_folder_hst+gal+'_'+case+'mas_f'+band+"w_wht.fits")
    except:
        hdul = fits.open(new_folder_hst+gal+'_'+case+'mas_f'+band+"lp_wht.fits")
    return hdul


# read sky value (precomputed)
def read_sky(gal, band):
    """ 'm3' is the same galaxy 3 (see readme) """
    if gal == "m3":
        gal = "3"
    elif gal == "m943":
        gal = "943"
    sky = np.loadtxt(new_folder_sky+"sky_"+gal+"_"+band+".txt")
    return sky


# gets the value of the central pixel in the segmentation mask
# pixels of the central object are identically labeled (except the border ones)
def get_key_from_segmentation_map(gal):
    segmentation_mask = fits.open(segmendir+'/udf_mosaic_'+gal+'_HST_SEGMAP.fits')
    N,M = segmentation_mask[1].data.shape
    key = segmentation_mask[1].data[N//2,M//2]
    return np.int(key)


# mask the region of the main object using the key on the central pixel
# the mask can also be configured to extend over two merging objects
def mask_main(gal, extra_key=None):
    """
    Reads the segmentation mask and returns a mask to separate the MAIN OBJECT
    ALTERNATIVELY, you can mask the TWO OBJECTS when there is a merger.

    For mergers use the parameter "extra_key", which you should set manually in
    the body of this function. You sould also define a custom name for the
    merger, for instance:
    gal = '3'   : masks only the main object based on the central pixel
    gal = 'm3'  : treats galaxy 3 as a merger of two objects

    Returns:
        a mask tagging pixels of the main object as TRUE (False elsewhere)
        the header of original segmentation mask (have info about coords.)

    NOTES:
        the mask returned is in the same grid of the segmentation maps
    """

    # Check if we want one these mergers encoded with the letter 'm' in front
    if gal == 'm3':
        gal = '3'
        extra_key = 24350 # identifies the second object in the segmentation map
    elif gal == 'm943':
        gal = '943'
        extra_key = 22951 # identifies the second object in the segmentation map

    # Read the segmentation mask
    segmentation_mask = fits.open(segmendir+'/udf_mosaic_'+gal+'_HST_SEGMAP.fits')
    seg_mask = segmentation_mask[1].data
    header = segmentation_mask[1].header

    # Get the key identifying the main object
    key = get_key_from_segmentation_map(gal)

    # Check if only one main object or two main objects (for mergers)
    if extra_key:
        index = (seg_mask == key) | (seg_mask == extra_key)
    else:
        index = (seg_mask == key)

    mask = seg_mask.copy()
    mask[index] = 1 # sky pixels (or other objects to be ignored)
    mask[~index] = 0 # source(s)
    mask = mask.astype(bool)
    return mask, header


# The boolean use of a single key misses some pixels in the edges of the object,
# so, the mask is enlarged by dilation in order to not loose any flux.
# Here we also return a version of the enlarged mask interpolated to the grid
# of the 160-band image to have all the data in the right frame of reference.
def enlarge_mask(gal, extra_key=None):
    band = '160'
    hdul = new_read_hst(gal, band, case='60')
    # this mask is in the grid of the segmentation maps
    mask, header = mask_main(gal, extra_key)
    # New creation of the mask including the dilation
    n = 3
    kernel = np.ones((n,n),np.uint8)
    mdil = cv2.dilate(mask.astype(float), kernel, iterations =1)
    # version of the mask mapped to the grid of the images
    new_mask = interpolate_mask(mdil, header, hdul)
    return new_mask


def interpolate_mask(mask, header, hdul):
    """
    Maps a mask from a native grid of pixels to a different one

    Receives the mask, its parent header, and the data over whose grid the
    mask is to be mapped, and uses interpolated-to-the-nearest.
    Returns a mask of 1's (and 0's) in the new grid of coordinates.
    """

    # here we do not need to oversample the data
    x,y = oversampled_positions(mask, oversampling=1)
    xx,yy = np.meshgrid(x, y)
    pixcrd = np.column_stack((np.ravel(yy), np.ravel(xx)))
    # we are converting the pixel coordinates to WCS
    world = get_wcs_coordinates(pixcrd, header)
    # go from WCS to pixel coordinates in the frame of ref of the other image
    pixcrd2 = get_pixel_coordinates(world, hdul[0].header)
    # interpolator is created in target frame of reference
    # new points are interpolated to the mask pixel they are nearest to (1 or 0)
    interp = interpolate.NearestNDInterpolator(pixcrd2, np.ravel(mask.T).T)
    # Let's apply the interpolating mask to the HST positions
    x,y = oversampled_positions(hdul[0].data, oversampling=1)
    xx,yy = np.meshgrid(x, y)
    pixcrd = np.column_stack((np.ravel(yy), np.ravel(xx)))
    new_mask = interp(pixcrd[:,0], pixcrd[:,1]).reshape(xx.shape).T
    new_mask = new_mask.astype(float)

    return new_mask


# (row, col) positions of a grid of oversampled pixels on top of a coarser grid
# 'oversample' gives the factor of oversampling
def oversampled_positions(image, oversampling):
    npix = image.shape[0]
    x = 0.5 + (1./oversampling)/2 +  np.arange(0, npix * oversampling, 1)/float(oversampling)
    npix = image.shape[1]
    y = 0.5 + (1./oversampling)/2 +  np.arange(0, npix * oversampling, 1)/float(oversampling)
    return x, y


# get WCS coordinates from pixel coordinates and header info
# WARNING: assuming FITS standard, pixels start at 1 (not at 0 !)
def get_wcs_coordinates(pixcrd, header):
    w = wcs.WCS(header)
    world = w.wcs_pix2world(pixcrd, 1)
    return world


# get pixel coordinates from WCS coordinates and header info
def get_pixel_coordinates(world,header):
    w = wcs.WCS(header)
    pixcrd = w.wcs_world2pix(world, 1)
    return pixcrd



