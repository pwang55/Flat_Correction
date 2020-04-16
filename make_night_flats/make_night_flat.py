"""

Usage:

    In script directory:
    $ python make_night_flat.py path_to_file/filelist.txt

    In files directory:
    $ python path_to_script/make_night_flat.py filelist.txt

This script take all the biasdarksub_xxx.fits in filelist.txt, cut out objects, normalized them and weighted-average all to make night_flat_asu1_feb.fits.


"""
import numpy as np
import sys
import os
from astropy.io import fits
from scipy.signal import savgol_filter
import re
import subprocess
from glob import glob

RN = 10.0
sig = 2.5

# if no argument given, return doc and exit
if len(sys.argv) != 2:
    print(__doc__)
    exit()

# absolute path of this script
script_path = sys.path[0]
# optic directory
optic_dir = os.environ['_90PRIME_OPTIC_DIR']



path_file = sys.argv[1]
# find input path and listfile name
filename = path_file.split('/')[-1]
if len(path_file.split('/')) > 1:
    filepath = path_file[:-len(filename)]
elif len(path_file.split('/')) == 1:
    filepath = ""
    
# read listfile
with open(path_file) as list0:
    frames0 = [ l.strip() for l in list0 ]

# figure out month
month = fits.getheader(filepath+frames0[0])['date'].split('-')[1]
# check to make sure all fits file have same month
months_check = []
for x in range(len(frames0)):
    mth = fits.getheader(filepath+frames0[x])['date'].split('-')[1]
    months_check.append(mth)
months_check = np.array(months_check)
if sum(months_check!=month) > 0:
    print('Some input files have different month, exiting!')
    exit()

if month == '02':
    month_dir = '/Feb/'
    month_str = 'feb'
elif month == '03':
    month_dir = '/Mar/'
    month_str = 'mar'


# For each frame, use sextractor to make a BACKGROUND_RMS as weight map, then use it to create segmentation map
for x in range(len(frames0)):

    base = frames0[x].replace('.fits','')
    if len(glob(filepath+base+'.bgrms.fits'))==0:   # if bgrms.fits doesn't exist
        subprocess.run(['sex' ,filepath+frames0[x], '-c', script_path+'/extra/make_bgrms_seg.config', '-PARAMETERS_NAME', script_path+'/extra/test.param', '-WEIGHT_TYPE', 'NONE', '-CHECKIMAGE_TYPE', 'BACKGROUND_RMS', '-CHECKIMAGE_NAME', filepath+base+'.bgrms.fits', '-FILTER_NAME', script_path+'/extra/gauss_4.0_7x7.conv', '-CATALOG_TYPE', 'NONE'])
    else:
        print(base+'.bgrms.fits already exists.')

    if len(glob(filepath+base+'.seg.fits'))==0:   # if seg.fits doesn't exist
        subprocess.run(['sex' ,filepath+frames0[x], '-c', script_path+'/extra/make_bgrms_seg.config', '-PARAMETERS_NAME', script_path+'/extra/test.param', '-WEIGHT_TYPE', 'MAP_RMS', '-WEIGHT_IMAGE', filepath+base+'.bgrms.fits', '-CHECKIMAGE_TYPE', 'SEGMENTATION', '-CHECKIMAGE_NAME', filepath+base+'.seg.fits', '-FILTER_NAME', script_path+'/extra/gauss_4.0_7x7.conv', '-CATALOG_TYPE', 'NONE'])
    else:
        print(base+'.seg.fits already exists.')


def combine_frames(filt='ASU1'):

    # Find frames that have the right filter
    frames = [ files for files in frames0 if (fits.getheader(filepath+files)['FILTER']==filt) ]
    print('Filter:\t'+filt)
    print('Creating Night Flat Using:')
    for x in range(len(frames)):
        print('\t',frames[x])

    # Find optic skymask file
    fskymask = optic_dir+month_dir+'optic_flat_'+filt.lower()+'_skymask.fits'
    skymask = fits.open(fskymask)
    print('Using SkyMask:\t'+'optic_flat_'+filt.lower()+'_skymask.fits')

    n = len(frames)
    with fits.open(filepath+frames[0]) as f:
        nx = f[1].header['naxis2']
        ny = f[1].header['naxis1']

    sky_max = np.zeros((n, 4), dtype=float)

    # Find normalization constant for each frame
    for i in range(n):
        temp_sky = []
        seg_file = frames[i].replace('.fits','.seg.fits')
        for j in range(1,17):
            with fits.open(filepath+frames[i]) as f:
                im = f[j].data
            with fits.open(filepath+seg_file) as sg:
                seg = sg[j].data
            im[seg>0] = 0.0             # Set all objects as zero
            ho = skymask[j].data > 0.0  # Find Skymask regions (the center flat part of filter)
            im[~ho] = 0.0               # Set outside sky region to zero so that we can find the median sky value as normalization
            msky = np.median(im[im>0])  # Find medium of the central sky as normalization constant
            temp_sky.append(msky)
        # For each filter, there are 4 sky value, pick the largest one as the single normalization constant for that fhilter
        for k in range(4):
            sky_max[i][k] = max(temp_sky[k*4], temp_sky[k*4+1], temp_sky[k*4+2], temp_sky[k*4+3])
        print('Sky of '+frames[i]+':'+'\t'+str(sky_max[i]))


    dat = np.zeros((n, nx, ny), dtype=float)
    dat2 = np.zeros((n, nx, ny), dtype=float)
    weights = np.zeros((n, nx, ny), dtype=float)

    hlist = []
    for j in range(1,17):   # iterate each chip
        for i in range(n):  # iterate each frame

            seg_file = frames[i].replace('.fits','.seg.fits')
            ho = skymask[j].data > 0.0

            with fits.open(filepath+seg_file) as sg:
                seg = sg[j].data

            print("Reading extension {} in {}".format(j, frames[i]))
            with fits.open(filepath+frames[i]) as d:
                im = d[j].data
                dat[i] = im/sky_max[i][(j-1)//4]
                h = seg > 0
                im[h] = 0.0
                im[~ho] = 0.0
                weights[i][True] = np.median(im[im>0])/np.std(im[im>0])
                dat[i][h] = 0.0
                weights[i][h] = 0.0
                dat2[i] = dat[i]*weights[i]

            std_j = np.std(dat, axis=0)
            med_j = np.median(dat, axis=0)

            # for each chip across different frames, if there is outlier, set it to zero
            for i in range(n):
                hs = (dat[i] < med_j + sig*std_j) & (dat[i] > med_j - sig*std_j)
                weights[i][~hs] = 0.0
                dat2[i][~hs] = 0.0

        hduI = fits.ImageHDU()
        datf = np.sum(dat2, axis=0)/np.sum(weights, axis=0)   # weighted average
        h = np.isnan(datf)
        datf[h] = 0.0
        hduI.data = datf
        hduI.header = fits.getheader(filepath+frames[0],j)
        hduI.header['BZERO'] = 0.
        hlist.append(hduI)

    hdu0 = fits.PrimaryHDU()
    hdu0.header = fits.getheader(filepath+frames[0])
    hlist.insert(0,hdu0)
    hduA = fits.HDUList(hlist)
    hduA.writeto(filepath+'night_flat_'+filt.lower()+'_'+month_str+'.fits')
    skymask.close()


combine_frames(filt='ASU1')
combine_frames(filt='ASU2')
combine_frames(filt='ASU3')
combine_frames(filt='ASU4')




