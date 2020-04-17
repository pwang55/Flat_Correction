"""
Usage:

    In script folder:
    $ python correct_flat_with_compare.py path_to_files/flatlist.txt

    In biassub files folder:
    $ python path_to_script/correct_flat_with_compare.py flatlist.txt


This script takes biassub_flat.fits and grab biassub_flat.compare.fits, use it to correct biassub_flat.fits,
output biassub_flat.corrected.fits

Input flatlist.txt should contain the original raw flat file names.


"""
import numpy as np
from astropy.io import fits
from glob import glob
from scipy.signal import savgol_filter as sf
from scipy.ndimage import median_filter as mf
import sys
import re
import os

median_size = 3

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
    filelists = [ filepath+'biassub_'+l.strip() for l in f ]

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

    print('Using Skymask:\t'+month_dir+'optic_flat_'+filt+'_skymask.fits')

    im = fits.open(filelists[i])
    compare = fits.open(comparelists[i])
    skymask = fits.open(fskymask)

    max_skys = []
    nx = im[1].data.shape[0]
    ny = im[1].data.shape[1]
    dat = np.zeros((16, nx, ny), dtype=float)
    for j in range(1,17):
        
        imj = im[j].data
        comparej = compare[j].data
        skymaskj = skymask[j].data

        mf_comparej = mf(comparej, size=median_size)
        dat[j-1] = imj/mf_comparej

        hsky = skymaskj > 0.0
        max_skys.append( np.median(dat[j-1][hsky]))

    max_sky2 = []
    for x in range(4):
        max_sky2.append( max(max_skys[x*4], max_skys[x*4+1], max_skys[x*4+2], max_skys[x*4+3]) )
    print(max_sky2)

    hlist = []
    for j in range(16):
        dat[j] = dat[j]/max_sky2[(j)//4]
        hduI = fits.ImageHDU()
        hduI.data = dat[j]
        hduI.header = im[j+1].header
        hlist.append(hduI)

    hdu0 = fits.PrimaryHDU()
    hdu0.header = im[0].header
    hlist.insert(0, hdu0)

    hduA = fits.HDUList(hlist)
    hduA.writeto(filelists[i].replace('.fits','.corrected.fits'))
    print('Corrected flat created:\t'+filelists[i].replace('.fits','.corrected.fits').split('/')[-1]+'\n')

    im.close()
    compare.close()
    skymask.close()


