from main import *

outdir = './data_example'
#galaxies = ["3","15","37","912","919","937","943","982","1002"]
galaxies = ["1","7","13"]


for gal in galaxies:

    for band in ["105","125","140","435","606","775","814","850"]:
        weight = read_hst_weight_maps(gal,band)

        newdata = signal_kernel_square(gal,band)
        filename = gal+'_'+band+'_kernel_squared.fits'
        fits.writeto(outdir+'/'+filename,newdata,weight[0].header,overwrite=True)



