"""
Usage:

    In script folder:
    $ python correct_flat_with_compare.py path_to_files/flatlist.txt

    In biassub files folder:
    $ python path_to_script/correct_flat_with_compare.py flatlist.txt


This script takes biassub_flat.fits and grab biassub_flat.compare.fits, use it to correct biassub_flat.fits,
output biassub_flat.corrected.fits

The file names in flatlist.txt should be the original raw flat file names.


"""
import numpy as np
from astropy.io import fits
from glob import glob
from scipy.signal import savgol_filter as sf
from scipy.ndimage import median_filter as mf
import sys
import re
import os
import multiprocessing as mp
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground, ModeEstimatorBackground, SExtractorBackground
from scipy.ndimage import gaussian_filter as gf

bkg_backsize = 20
bkg_fsize = 15
bkg_estimator = MedianBackground()
sigma_clip = SigmaClip(sigma = 3.)
smoothing_method = 'gaussian'

# Smooth background after creating background
median_size = 3
gaussian_size = 3


if len(sys.argv)!=2:
    print(__doc__)
    exit()


path_file = sys.argv[1]
filename = path_file.split('/')[-1]
if len(path_file.split('/'))>1:
    filepath = path_file[:-len(filename)]
else:
    filepath = ""

with open(path_file) as f:
    filelists = [filepath + 'biassub_' + l.strip() for l in f]
    # filelists = [filepath + l.strip() for l in f]

comparelists = [ l.replace('.fits','.compare.fits') for l in filelists ]


n = len(filelists)

script_path = sys.path[0]
optic_dir = os.environ['_90PRIME_OPTIC_DIR']
bias_dark_flats_dir = os.environ['_90PRIME_BIAS_DARK_FLATS_DIR']


for i in range(n):
    print('Correct\t'+filelists[i].split('/')[-1]+' with\t'+comparelists[i].split('/')[-1])

    filt = fits.getheader(filelists[i])['filter'].lower()
    month = fits.getheader(filelists[i])['date'].split('-')[1]
    if month == '02':
        month_dir = '/Feb/'
        month_str = 'feb'
    elif month == '03':
        month_dir = '/Mar/'
        month_str = 'mar'
    fskymask = optic_dir+month_dir+'optic_flat_'+filt+'_skymask.fits'

    im = fits.open(filelists[i])
    compare = fits.open(comparelists[i])

    nx = im[1].data.shape[0]
    ny = im[1].data.shape[1]
    dat = np.zeros((16, nx, ny), dtype=float)


    hlist = []
    for j in range(1, 17):
                
        imj = im[j].data
        comparej = compare[j].data
        # skymaskj = skymask[j].data

        compare_bkg = Background2D(comparej, (bkg_backsize, bkg_backsize), filter_size=(bkg_fsize, bkg_fsize), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

        if smoothing_method == 'gaussian':
            compare_final = gf(compare_bkg.background, gaussian_size)
        elif smoothing_method == 'median':
            compare_final = mf(compare_bkg.background, median_size)

        icorrected = imj / compare_final

        hduI = fits.ImageHDU()
        hduI.data = icorrected
        hduI.header = im[j].header
        hlist.append(hduI)

    hdu0 = fits.PrimaryHDU()
    hdu0.header = im[0].header
    hlist.insert(0, hdu0)

    hduA = fits.HDUList(hlist)
    hduA.writeto(filelists[i].replace('.fits','.corrected.fits'), overwrite=True)
    print('Corrected flat created:\t'+filelists[i].replace('.fits','.corrected.fits').split('/')[-1]+'\n')

    im.close()
    compare.close()



