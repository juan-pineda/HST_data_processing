from main import *

outdir = './results_example'
galaxies = ["3","15","37","912","919","937","943","982","1002","m3","m943"]

# NOTE:
# Galaxies (3, 943) are mergers. If we want to retain the main object only we
# use the number of the galaxy (e.g. 3). We use an extra 'm' in front (e.g. m3)
# if we want to model the flux of the merging objects as a single thing
# This is very relevant for the segmentation masks.

for gal in galaxies:
    for band in ["105","125","140","435","606","775","814","850"]:
        weight = read_hst_weight_maps(gal,band)

        # low resolution version of the image
        newdata = signal_low_resolution(gal,band)
        filename = gal+'_'+band+'_lowres.fits'
        fits.writeto(outdir+'/'+filename,newdata,weight[0].header,overwrite=True)

        # low resolution version of the variance map
        newdata = variance_low_resolution(gal,band)
        filename = gal+'_'+band+'_variance_lowres.fits'
        fits.writeto(outdir+'/'+filename,newdata,weight[0].header,overwrite=True)




