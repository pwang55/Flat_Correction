"""
Usage:

    In script folder:
    $ python compare_biassub_flat2nights.py path_to_file/flatlist.txt

    In biassub flats folder:
    $ python path_to_script/compare_biassub_flat2nights.py flatlist.txt

This script takes biassub_flat.fits, compare it to night flats, and generate a .compare.fits that is the ratio of 
biassub_flat.fits to night flats, used later for correcting twilight flats.

The input flatlist.txt should be the original one containing raw flat file names, not the biassub ones or compare/corrected ones.


"""
import numpy as np
from astropy.io import fits
import sys
from glob import glob
import os


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
    filelists = [ filepath+'biassub_'+l.strip() for l in f ]

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

    ffiltermask = optic_dir+month_dir+'optic_flat_'+filt+'_filtermask.fits'
    fskymask = optic_dir+month_dir+'optic_flat_'+filt+'_skymask.fits'
    fnightflat = bias_dark_flats_dir+'/nightflats'+month_dir+'night_flat_'+filt+'.fits'

    print('\nFiltermask:\t'+'optic_flat_'+filt+'_filtermask.fits')
    print('Skymask:\t'+'optic_flat_'+filt+'_skymask.fits')
    print('Night flat:\t'+month_dir+'night_flat_'+filt+'.fits\n')

    filtermask = fits.open(ffiltermask)
    skymask = fits.open(fskymask)
    nightflat = fits.open(fnightflat)

    hlist = []
    sky_count = []

    for j in range(1,17):   # iterate each chip
        nightj = nightflat[j].data
        skymaskj = skymask[j].data
        filtermaskj = filtermask[j].data
        imj = im0[j].data
        
        dat0 = np.zeros((nightj.shape[0], nightj.shape[1]), dtype='float32')
        h0 = filtermaskj > 0.0  # filtermask > 0
        dat0[h0] = imj[h0]/nightj[h0]
        msky = np.median(dat0[skymaskj>0])
        sky_count.append(msky)

    sky_2 = []
    for k in range(4):
        max_sky2 = max(sky_count[k*4], sky_count[k*4+1], sky_count[k*4+2], sky_count[k*4+3])
        sky_2.append(max_sky2)

    print(sky_2)

    for j in range(1,17):
        nightj = nightflat[j].data
        skymaskj = skymask[j].data
        filtermaskj = filtermask[j].data
        imj = im0[j].data 

        dat = np.zeros((nightj.shape[0], nightj.shape[1]), dtype='float32')
#        h = filtermaskj > 0.0
#        dat[h] = imj[h]/nightj[h]/sky_2[(j-1)//4]
        dat = imj/nightj/sky_2[(j-1)//4]

        hduI = fits.ImageHDU()
        hduI.data = dat
        hduI.header = im0[j].header
        hlist.append(hduI)

    hdu0 = fits.PrimaryHDU()
    hdu0.header = im0[0].header
    hlist.insert(0,hdu0)

    hduA = fits.HDUList(hlist)
    hduA.writeto(filelists[i].replace('.fits','.compare.fits'))
    print('Compare file created:\t'+filelists[i].replace('.fits','.compare.fits').split('/')[-1])

    im0.close()
    filtermask.close()
    skymask.close()
    nightflat.close()



