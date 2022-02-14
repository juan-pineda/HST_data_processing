# HST_data_processing

This module deals with HST data, basically images and weight maps defined as
the inverse of the variance in the flux. Segmentation masks are also needed as
input for the processing.   

The purpose is to create for a given galaxy a set of standard maps across all
bands. This is done by functions in the module "main.py", which basically make
use of the segmentation masks to retain only the main object and lower its 
resolution to the poorest among the different bands (160 microns).   

The module "functions.py" is in charge of lower-level tasks. On it you must
define the paths to the folders with the HST data, measurements of the sky
in the different images, and the folder with the segmentation maps. For tests
we have provided here some data for a subsample of galaxies, which should allow
you to run the "test.py" file.  

### Format of HST filenames is, for instance for galaxy 872:
- "872_60mas_f814w_sci.fits" (science image)  
- "872_60mas_f814w_wht.fits" (weight image)  
- "872_60mas_f850lp_sci.fits" (only this band has suffix lp)
- There are also 30mas images

### Segmentation masks
These were provided by Thierry Contini. They do not have the same orientation
as the images, which we solve by an interpolation. These interpolation routines
are quite powerful, as they allow the mapping between data coming from different
instruments, as we explain somewhere else for matching HST images on top of OII
line fluxes from MUSE.  

In the masks different objects are represented by different
numbers. They are centered in the main object, so the central value gives the
correct mapping to retain only this primary object. Nevertheless, for a weird
reason the values of the pixels in the edge of the main object change, so one
may loose some flux. To avoid it we enlarge the masks by performing a dilation.  


### On SKY and PSF resolution
The SKY has to be measured in advance. This is explained somewhere else, here we
directly provide those measurements as text files in folder NEW_SKY for the data
we were using. This was measured for the 60mas images.  

The spatial resolution was measured in advance as well for the 60mas images.
This is done by fitting gaussian profiles to 5 stars in each band, and will be
explained in more detail somewhere else. The estimated resolutions are manually
copied inside the file "main.py".  



