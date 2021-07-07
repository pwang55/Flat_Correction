"""
Usage:

    In script folder:
    $ python flat_sub_bias.py path_to_files/flatlist.txt

    In flat files folder:
    $ python path_to_script/flat_sub_bias.py flatlist.txt

This script takes flat files and subtract it with bias.


"""
import numpy as np
from astropy.io import fits
import os
from scipy.signal import savgol_filter
import re
import sys
from glob import glob


apply_overscan = True

if len(sys.argv) !=2:
    print(__doc__)
    exit()


# Get absolute path of this script
script_path = sys.path[0]
# Get optics path and bias, dark, flat dir
optic_dir = os.environ['_90PRIME_OPTIC_DIR']
bias_dark_flats_dir = os.environ['_90PRIME_BIAS_DARK_FLATS_DIR']

path_file = sys.argv[1]
filename = path_file.split('/')[-1]
if len(path_file.split('/'))>1:
    filepath = path_file[:-len(filename)]
else:
    filepath = ""

# Read in filenames and store them in a list with full path to the files
with open(path_file) as f:
    filelists = [ filepath+l.strip() for l in f ]


# Define function to find overscan region; these n number will be the exact pixel number (start 0) of where biassec start and finish, so if using these in argument, the end has to add 1
# 90prime data NAXIS1 = y, NAXIS2 = x, so in header [NAXIS1, NAXIS2], but when using loaded array, dat[naxis2, naxis1]
# default to use the end of naxis1 = y as overscan
def find_overscan(biassec, datasec):
    biasstr = re.split('[\[:,\]]', biassec)
    datastr = re.split('[\[:,\]]', datasec)
    nyb0 = int(biasstr[1])-1
    nyb1 = int(biasstr[2])-1
    nxb0 = int(biasstr[3])-1
    nxb1 = int(biasstr[4])-1
    nyd0 = int(datastr[1])-1
    nyd1 = int(datastr[2])-1
    nxd0 = int(datastr[3])-1
    nxd1 = int(datastr[4])-1
    return nxb0, nxb1, nyb0, nyb1, nxd0, nxd1, nyd0, nyd1


n = len(filelists)

# Get necessary dimensions

with fits.open(filelists[0]) as f:
    namp = len(f)                    # number of amplifiers+1 (including Primary HDU)
    ovsy = int(f[1].header['ovrscan1'])
    ovsx = int(f[1].header['ovrscan2'])
    ny = f[1].header['naxis1'] - ovsy
    nx = f[1].header['naxis2'] - ovsx
    month = f[0].header['date'].split('-')[1]
    if month == '02':
        month_string = 'feb'
    elif month == '03':
        month_string = 'mar'

if apply_overscan == True:
    biasfile = bias_dark_flats_dir + '/bias/bias.fits'
else:
    biasfile = bias_dark_flats_dir + '/bias/bias_' + month_string + '.fits'


for x in range(n):

    print('Subtract bias from:\t'+filelists[x].split('/')[-1])
    hlist = []
    with fits.open(filelists[x]) as im0:
        for j in range(1, namp):
            nxb0, nxb1, nyb0, nyb1, nxd0, nxd1, nyd0, nyd1 = find_overscan(im0[j].header['biassec'], im0[j].header['datasec'])
            im = im0[j].data[nxd0:nxd1+1, nyd0:nyd1+1]

            if apply_overscan == True:
                overscan0 = np.average(im0[j].data[nxd0:nxd1+1, nyb0:nyb1+1], axis=1)
                overscan1 = savgol_filter(overscan0, 201, 2)
                im = im - overscan1[:, None] - fits.getdata(biasfile, j)
            else:
                im = im - fits.getdata(biasfile, j)

            dat = im

            hduI = fits.ImageHDU()
            hduI.data = dat
            hduI.header = im0[j].header
            hduI.header['BZERO'] = 0.
            hlist.append(hduI)

        hdu0 = fits.PrimaryHDU()
        hdu0.header = im0[0].header
        hlist.insert(0,hdu0)
    hduA = fits.HDUList(hlist)
    hduA.writeto(filepath+'biassub_'+filelists[x].split('/')[-1])




