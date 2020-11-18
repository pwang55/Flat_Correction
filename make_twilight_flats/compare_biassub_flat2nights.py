"""
Usage:

    In script folder:
    $ python compare_biassub_flat2nights.py path_to_file/flatlist.txt

    In biassub flats folder:
    $ python path_to_script/compare_biassub_flat2nights.py flatlist.txt

This script takes biassub_flat.fits, compare it to night flats, and generate a .compare.fits that is the ratio of 
biassub_flat.fits to night flats, used later for correcting twilight flats.

The file names in flatlist.txt should be the original raw flat file names.


"""
import numpy as np
from astropy.io import fits
import sys
from glob import glob
import os
from scipy.ndimage import median_filter as mf
import time

median_fulter_size = 64


t1 = time.time()

if len(sys.argv) !=2:
    print(__doc__)
    exit()


path_file = sys.argv[1]
filename = path_file.split('/')[-1]
if len(path_file.split('/'))>1:
    filepath = path_file[:-len(filename)]
else:
    filepath = ""

# Get absolute path of this script
script_path = sys.path[0]
# Get optic directory and bias dark flats directory
optic_dir = os.environ['_90PRIME_OPTIC_DIR']
bias_dark_flats_dir = os.environ['_90PRIME_BIAS_DARK_FLATS_DIR']

# read filename into a list and filename contains full path to it
with open(path_file) as f:
    filelists = [filepath + 'biassub_' + l.strip() for l in f]
    # filelists = [filepath + l.strip() for l in f]
    


n = len(filelists)

for i in range(n):
    print('Making Compare file for:\t'+filelists[i].split('/')[-1])
    im0 = fits.open(filelists[i])
    filt = im0[0].header['filter'].lower()
    month = im0[0].header['date'].split('-')[1]
    if month == '02':
        month_dir = '/Feb/'
        month_str = 'feb'
    elif month == '03':
        month_dir = '/Mar/'
        month_str = 'mar'


    fnightflat = bias_dark_flats_dir + '/nightflats' + month_dir + 'night_flat_' + filt + '_' + month_str + '.fits'
    print('Night flat:\t'+month_dir+'night_flat_'+filt+'.fits\n')
    nightflat = fits.open(fnightflat)

    hlist = []
    sky_count = []



    for j in range(1,17):
        nightj = nightflat[j].data
        imj = im0[j].data 

        dat = np.zeros((nightj.shape[0], nightj.shape[1]), dtype='float32')
        dat = imj/nightj
        # mdat = mf(dat, median_fulter_size)
        mdat = dat
        hduI = fits.ImageHDU()
        hduI.data = mdat
        hduI.header = im0[j].header
        hlist.append(hduI)

    hdu0 = fits.PrimaryHDU()
    hdu0.header = im0[0].header
    hlist.insert(0,hdu0)

    hduA = fits.HDUList(hlist)
    hduA.writeto(filelists[i].replace('.fits','.compare.fits'), overwrite=True)
    print('Compare file created:\t'+filelists[i].replace('.fits','.compare.fits').split('/')[-1])

    im0.close()
    nightflat.close()
    t2 = time.time()
    print('Time: {}\n'.format(t2 - t1))


