**Authors:**  
- Juan Carlos Basto Pineda (juan.basto.pineda@gmail.com)
- Benoit Epinat

**Updated:** 04/10/2022

# HST\_data\_processing

This pipeline creates homogeneous data products for a given field across
different bands of HST. Those data products can be images or weight maps, the
latter defined as the inverse of the variance of the flux.

By homogeneous we mean:

- Data is smoothed to the same resolution in all bands
- That resolution is the one of the poorest (160 band)
- All maps are defined exactly in the very same pixels
- Those pixels are given by the segmentation maps

To quickly grasp how this repository works you can check now the example.py file
and come back here after that.  

---

This is a whole repo on its own, in the sense that serves to a specific purpose
that it is useful in itself. Yet it is important to be aware that there is a
[second repository](https://github.com/juan-pineda/Models_HST_MUSE) coming
which recycles a number of the functions given here. We keep them separated for
several reasons:

- In the second one we have a very specific goal: to create approximate models
of the OII flux of a given object by combining two HST images.
- There we do not work with weight maps
- We will consider the possibility of using images with different platescales,
i.e., 30mas and 60mas, while here we deal with 60mas images only

Overall, if your work is more related to the second repo, we suggest you start
here anyways, as this pipeline is simpler and will give you an ideal background
to move on to the other project, whose structure is more challenging to follow
from scratch.

# INPUT DATA

For tests we have provided some data for a subsample of galaxies, which should
allow you to run the example.py file, after making sure to correctly set the
paths given inside functions.py.


## new\_stamps

Contains images and weight maps in 9 bands for different HST fields. Here we
only deal with 60mas data.  

Format of HST filenames, for instance for galaxy 37 and band 814, is:
- 37\_60mas\_f814w\_sci.fits (science image)
- 37\_60mas\_f814w\_wht.fits (weight image)
- 37\_60mas\_f850lp\_sci.fits (only this band has suffix lp)

The spatial resolution was measured in advance for the 60mas images in all bands.
This is done by fitting gaussian profiles to 5 stars in each band, as explained
in the [second repository](https://github.com/juan-pineda/Models_HST_MUSE)
The estimated resolutions are manually declared inside the file main.py.

## segment\_maps

Segmentation maps needed to separate the object from the background. Note that
these maps do not have the same orientation as the images. We fix this by
incorporating some non-trivial interpolation routines that allow the mapping 
of data coming from different instruments to a common grid.  

In the segmentation maps different objects are represented by different numbers.
They are centered in the main object, so, the central value gives the correct
key to retain only this primary object. Nevertheless, for an unknown reason 
pixels in the edge of an object exhibit different values. To avoid missing those
pixels we create initial masks using the central key and then enlarge them by
dilation.  

In some cases a field is targeting a merging system and each galaxy will be
identified by a distinct key. As the authomatic mask extraction will only retain
the pixels from the primary object, in case you want to create the maps for both
objects together a manual step must be followed, to tell the code what is the
key of the second object. Instructions are given in the file functions.py

## NEW\_SKY

For each galaxy-band combination there is a text file with the estimated value
of the sky and a plot of the histogram of pixels in the background from which
the sky was computed. This files are precomputed, the way to do it is explained
in the [second repository](https://github.com/juan-pineda/Models_HST_MUSE) 


# Example
The test case in example.py. has been checked and must run smoothly if you have
the necessary libraries installed and a healthy python distribution. Make sure
you can run it, you understand what it did, and explore the results to see if
they look as expected before building on top of this example or expand its use
to a larger sample of galaxies.


