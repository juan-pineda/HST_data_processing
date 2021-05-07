
# This script is based in the notebook called `calc_sky_using_segmentation_mask.ipynb`

import sys
sys.path.append("../")
from functions_segmentation_masks import *


#galaxies = ["3","15","37","912","919","937","943","982","1002"]
galaxies = ["1","7","13"]
bands = ["105","125","140","160","435","606","775","814","850"]

#galaxies = [str(gal) for gal in [1,2,4,5,7,8,11,12,13,14]]

for gal in galaxies:
    print("working in gal ",gal)
    hdul1, hdul2 = read_muse(gal)
    for band in bands:
        hdul3 = new_read_hst(gal,band,case='60')
        mask, header, seg_mask = mask_sky(gal)
        new_mask = interpolate_mask(mask,header,hdul3)

        if band in ["435","606","775","814","850"]:
            sky, std = get_sky_n_fluctuations(hdul3,new_mask,graph=True,filename="histogram_"+gal+"_"+band+".png")
        elif band in ["105","125","140","160"]:
            erosion = erode_sky(new_mask,n=5,iterations=5)
            sky, std = get_sky_n_fluctuations(hdul3,erosion,graph=True,filename="histogram_"+gal+"_"+band+".png")
        filename = "sky_"+gal+"_"+band+".txt"
        np.savetxt(filename,np.array([sky]))

        hdul3.close()
    hdul1.close()
    hdul2.close()


