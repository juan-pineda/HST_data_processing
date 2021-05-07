from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy import wcs
import scipy.interpolate
from numpy.fft import fftshift, rfft2, irfft2
from shutil import copyfile
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
import os
import sys
import cv2
from scipy import interpolate

# new_folder_hst must point to the HST data
# format is:
# 872_60mas_f814w_sci.fits
# 872_60mas_f814w_wht.fits
# 872_60mas_f850lp_sci.fits
new_folder_hst = '/media/juan/Pesquisa/DATOS_MUSE/new_stamps/'

# Folder with the measurements of the sky in the background of each image
new_folder_sky = "./NEW_SKY/"

# Here are the segmentation masks
segmendir = '/home/juan/PROYECTOS/MUSE/Data_from_Contini/mosaics'



# Read HST image
def new_read_hst(gal, band, case='60'):
    """
    'case' stands for pixel size in mas, must be: '30' or '60'
    """
    try:
         hdul = fits.open(new_folder_hst+gal+'_'+case+'mas_f'+band+"w_sci.fits")
    except:
        hdul = fits.open(new_folder_hst+gal+'_'+case+'mas_f'+band+"lp_sci.fits")
    return hdul


# Read weight map
def read_hst_weight_maps(gal, band, case='60'):
    """
    'case' stands for pixel size in mas, must be: '30' or '60'
    """
    try:
         hdul = fits.open(new_folder_hst+gal+'_'+case+'mas_f'+band+"w_wht.fits")
    except:
        hdul = fits.open(new_folder_hst+gal+'_'+case+'mas_f'+band+"lp_wht.fits")
    return hdul


# Read sky value
def read_sky(gal, band):
    sky = np.loadtxt(new_folder_sky+"sky_"+gal+"_"+band+".txt")
    return sky


# It gets the central pixel of the segmentation mask, marking the main object
def get_key_from_segmentation_map(gal):
    segmentation_mask = fits.open(segmendir+'/udf_mosaic_'+gal+'_HST_SEGMAP.fits')
    N,M = segmentation_mask[1].data.shape
    key = segmentation_mask[1].data[N//2,M//2]
    return np.int(key)


# Mask the region of the main object using the key on the central pixel
# extra key is intended for mergers, when you want to retain two objects
# gal = 'm3': treats galaxy 3 as a merger of two objects
# gal = 'm943': treats galaxy 943 as a merger of two objects
def mask_main(gal,extra_key=None):
    """
    Reads the segmentation mask, and returns a mask to separate the MAIN OBJECT
    Returns:
    mask: mask indicating the main object with TRUEs (and False elsewhere)
    header: header from the original segmentation mask (have info about coords.)
    """
    # Check if we want one these mergers, encoded with the letter 'm' in front
    if gal == 'm3':
        gal = '3'
        extra_key = 24350
    elif gal == 'm943':
        gal = '943'
        extra_key = 22951

    # Read the segmentation mask
    segmentation_mask = fits.open(segmendir+'/udf_mosaic_'+gal+'_HST_SEGMAP.fits')
    seg_mask = segmentation_mask[1].data
    header = segmentation_mask[1].header

    # Get the key identifying the main object
    key = get_key_from_segmentation_map(gal)

    # Check for main object (or two main objects)
    if extra_key:
        index = (seg_mask == key) | (seg_mask == extra_key)
    else:
        index = (seg_mask == key)

    mask = seg_mask.copy()
    mask[index] = 1 # sky pixels
    mask[~index] = 0 # sources
    mask = mask.astype(bool)
    return mask, header


# As boolean search for single key misses some pixels in the edges,
# the maks is enlarged by dilation
def enlarge_mask(gal,extra_key=None):
    band = '160'
    hdul = new_read_hst(gal,band,case='60')
    mask, header = mask_main(gal,extra_key)
    # New creation of the mask including the dilation
    n = 3
    kernel = np.ones((n,n),np.uint8)
    mdil = cv2.dilate(mask.astype(float),kernel,iterations =1)
    new_mask = interpolate_mask(mdil,header,hdul)
    return new_mask


def interpolate_mask(mask,header,hdul):
    """
    Receives the mask, its parent header, and the image
    over which grid the mask is to mapped;
    mask is interpolated to the nearest in the frame of
    reference of the image
    Returns a mask of 1's (and 0's) in these new coordinates
    """
    x,y = oversampled_positions(mask, oversampling=1)
    xx,yy = np.meshgrid(x,y)
    pixcrd = np.column_stack((np.ravel(yy),np.ravel(xx)))
    world = get_wcs_coordinates(pixcrd, header)
    # interpolator is created in target frame of reference
    pixcrd2 = get_pixel_coordinates(world, hdul[0].header)
    interp = interpolate.NearestNDInterpolator(pixcrd2, np.ravel(mask.T).T)
    # Let's apply the mask to the HST image
    x,y = oversampled_positions(hdul[0].data, oversampling=1)
    xx,yy = np.meshgrid(x,y)
    pixcrd = np.column_stack((np.ravel(yy),np.ravel(xx)))
    new_mask = interp(pixcrd[:,0],pixcrd[:,1]).reshape(xx.shape).T
    new_mask = new_mask.astype(float)
    return new_mask


# row, col positions of a grid of pixels oversampled on top of a coarser grid
# oversample gives the factor of oversampling
def oversampled_positions(image, oversampling):
    npix = image.shape[0]
    x = 0.5 + (1./oversampling)/2 +  np.arange(0, npix * oversampling, 1)/float(oversampling)
    npix = image.shape[1]
    y = 0.5 + (1./oversampling)/2 +  np.arange(0, npix * oversampling, 1)/float(oversampling)
    ###
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



