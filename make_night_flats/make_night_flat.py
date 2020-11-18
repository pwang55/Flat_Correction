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
from astropy.io import fits, ascii
# from scipy.signal import savgol_filter
from scipy.ndimage import median_filter as mf
from scipy.ndimage import gaussian_filter as gf
import re
import subprocess
from glob import glob
import multiprocessing as mp
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground, ModeEstimatorBackground, SExtractorBackground
# from astropy.convolution import Gaussian2DKernel
# from astropy.stats import gaussian_fwhm_to_sigma
# from photutils import detect_sources
# from photutils import detect_threshold
import collections
from astropy.table import Table

backsize = 64
# threshold_sig = 2.0
sigma_clip = SigmaClip(sigma = 3.)  # sigmaclip for background estimation
bkg_filter_size = 7
bkg_estimator = MedianBackground()
# bkg_estimator = ModeEstimatorBackground()
# sigma = 3.0 * gaussian_fwhm_to_sigma    # FWHM of kernel
# kernel = Gaussian2DKernel(sigma, x_size = 3, y_size = 3)
# kernel.normalize()


sig = 1.0   # sigma for weight average rejection
gf_size = 5
seg_fillvalue = 50


# Mask for chip j=9
xi = [1632, 1631, 1630, 1630, 1632, 1633]
xf = [1638, 1639, 1639, 1638, 1638, 1637]
yi = [1289, 1291, 1294, 1295, 1300, 1305]
yf = [1291, 1294, 1295, 1300, 1305, 2048]
# Set any additional width in pixels to extend from the above mask
xmargin = 2
ymargin = 1

# x, y grid for circle masking bright stars' outer light
xgrid = np.linspace(0, 2015, 2016, dtype = int)
ygrid = np.linspace(0, 2047, 2048, dtype = int)
xgrid, ygrid = np.meshgrid(xgrid, ygrid)



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

    if len(glob(filepath + base + '.seg.fits')) == 0:
        subprocess.run(['sex', filepath + frames0[x], '-c', script_path + '/extra/make_bgrms_seg.config', '-PARAMETERS_NAME', script_path + '/extra/test.param', \
        '-WEIGHT_TYPE', 'BACKGROUND', '-CHECKIMAGE_TYPE', 'SEGMENTATION', '-CHECKIMAGE_NAME', filepath + base + '.seg.fits', \
        '-FILTER_NAME', script_path + '/extra/gauss_3.0_7x7.conv', '-CATALOG_TYPE', 'NONE'])
    else:
        print(base + '.seg.fits exists')

# Function for running without multiprocessing
def find_normalization(j, f, sg, sky):

    im = f
    seg = sg

    # set the long straignt line to zero in seg
    if j == 9:
        for m in range(len(xi)):
            seg[yi[m] - ymargin:yf[m] + ymargin, xi[m] - xmargin:xf[m] + xmargin] = 0
    
    seg[seg > 0] = seg_fillvalue
    seg = gf(seg, gf_size)
    skymask = sky

    im[(seg > 0) | (skymask == 0)] = 0.0  # Set all objects as zeros
    msky = np.median(im[im > 0])
    return msky

# Function that run 4 chips for multiprocessing
def find_normalization_4(ji, jf, framei, sky, child):

    temp_sky1 = []
    ims = fits.open(framei)
    sgs = fits.open(framei.replace('.fits', '.seg.fits'))

    skys = fits.open(sky)
    for j in range(ji, jf + 1):
        im = ims[j].data
        seg = sgs[j].data

        # set the long straignt line to zero in seg
        if j == 9:
            for m in range(len(xi)):
                seg[yi[m] - ymargin:yf[m] + ymargin, xi[m] - xmargin:xf[m] + xmargin] = 0

        seg[seg > 0] = seg_fillvalue
        seg = gf(seg, gf_size)
        skymask = skys[j].data

        # bkg = Background2D(im, (backsize, backsize), filter_size=(bkg_filter_size, bkg_filter_size), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        # ibkg = bkg.background
        im[(seg > 0) | (skymask == 0)] = 0.0  # Set all objects as zeros
        # ibkg[skymask == 0] = 0.0
        # im = mf(im, size=mf_size)
        temp_sky1.append(np.median(im[im > 0]))  # Find medium of the central sky as normalization constant
        # temp_sky1.append(np.median(ibkg[ibkg > 0]))
    # return temp_sky1
    ims.close()
    sgs.close()
    skys.close()
    child.send(temp_sky1)
    child.close()

 
def combine_nightflat(j, n, nx, ny, sky, frames, sky_max):

    dat = np.zeros((n, nx, ny), dtype=float)
    dat2 = np.zeros((n, nx, ny), dtype=float)
    weights = np.zeros((n, nx, ny), dtype=float)

    for i in range(n):  # iterate each frame

        seg_file = frames[i].replace('.fits','.seg.fits')
        ho = sky > 0.0

        with fits.open(filepath + seg_file) as sg:
            seg = sg[j].data

        # set the long straignt line to zero in seg
        if j == 9:
            for m in range(len(xi)):
                seg[yi[m] - ymargin:yf[m] + ymargin, xi[m] - xmargin:xf[m] + xmargin] = 0

        seg[seg > 0] = seg_fillvalue
        seg = gf(seg, gf_size)

        print("Reading extension {} in {}".format(j, frames[i]))
        with fits.open(filepath + frames[i]) as d:
            im = d[j].data
            dat[i] = im / sky_max[i][(j - 1) // 4]
            h = seg > 0
            im[h] = 0.0
            im[~ho] = 0.0
            # weights[i][True] = np.median(im[im>0])/np.std(im[im>0])
            weights[i][True] = 1.0 / np.std(im[im > 0])** 2
            dat[i][h] = np.nan
            weights[i][h] = np.nan
            dat2[i] = dat[i] * weights[i]

    std_j = np.nanstd(dat, axis=0)
    med_j = np.nanmedian(dat, axis=0)
    # std_j = np.std(dat, axis=0)
    # med_j = np.median(dat, axis=0)

    # for each chip across different frames, if there is outlier, set it to zero
    for i in range(n):
        hs = (dat[i] < med_j + sig * std_j) & (dat[i] > med_j - sig * std_j)
        weights[i][~hs] = 0.0
        dat2[i][~hs] = 0.0

    hduI = fits.ImageHDU()
    datf = np.sum(dat2, axis=0) / np.sum(weights, axis=0)  # weighted average
    # datf = np.nanmedian(dat, axis=0)    # median
    h = np.isnan(datf)
    datf[h] = 0.0
    hduI.data = datf
    hduI.header = fits.getheader(filepath + frames[0], j)
    hduI.header['BZERO'] = 0.
    return hduI
    # hlisti.append(hduI)



def combine_nightflat_4(ji, jf, n, nx, ny, fskymask, frames, sky_max, child):

    dat = np.zeros((n, nx, ny), dtype=float)
    dat2 = np.zeros((n, nx, ny), dtype=float)
    weights = np.zeros((n, nx, ny), dtype=float)

    hlisti = []
    skymask = fits.open(fskymask)

    for j in range(ji, jf + 1):
        sky = skymask[j].data
        for i in range(n):  # iterate each frame

            seg_file = frames[i].replace('.fits','.seg.fits')
            ho = sky > 0.0

            with fits.open(filepath + seg_file) as sg:
                seg = sg[j].data

            # set the long straignt line to zero in seg
            if j == 9:
                for m in range(len(xi)):
                    seg[yi[m] - ymargin:yf[m] + ymargin, xi[m] - xmargin:xf[m] + xmargin] = 0

            seg[seg > 0] = seg_fillvalue
            seg = gf(seg, gf_size)

            print("Reading extension {} in {}".format(j, frames[i]))
            with fits.open(filepath + frames[i]) as d:
                im = d[j].data
                # bkg = Background2D(im, (backsize, backsize), filter_size=(bkg_filter_size, bkg_filter_size), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
                dat[i] = im / sky_max[i][(j - 1) // 4]
                h = seg > 0
                im[h] = 0.0
                im[~ho] = 0.0
                # weights[i][True] = np.median(im[im>0])/np.std(im[im>0])
                weights[i][True] = 1.0 / np.std(im[im > 0])** 2
                dat[i][h] = np.nan
                weights[i][h] = np.nan
                dat2[i] = dat[i] * weights[i]

        std_j = np.nanstd(dat, axis=0)
        med_j = np.nanmedian(dat, axis=0)
        # std_j = np.std(dat, axis=0)
        # med_j = np.median(dat, axis=0)

        # for each chip across different frames, if there is outlier, set it to zero
        for i in range(n):
            hs = (dat[i] < med_j + sig * std_j) & (dat[i] > med_j - sig * std_j)
            weights[i][~hs] = 0.0
            dat2[i][~hs] = 0.0

        hduI = fits.ImageHDU()
        datf = np.sum(dat2, axis=0) / np.sum(weights, axis=0)  # weighted average
        # datf = np.nanmedian(dat, axis=0)    # median
        h = np.isnan(datf)
        datf[h] = 0.0
        hduI.data = datf
        hduI.header = fits.getheader(filepath + frames[0], j)
        hduI.header['BZERO'] = 0.
        # return hduI
        hlisti.append(hduI)
    child.send(hlisti)
    child.close()




# ========================================================================================================================

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
        f = fits.open(filepath + frames[i])
        sg = fits.open(filepath + seg_file)

        # temp_sky = [find_normalization(x, f[x].data, sg[x].data, skymask[x].data) for x in range(1,17)]
        # temp_sky1 = find_normalization_4(1, 4, frames[i], fskymask)
        # temp_sky2 = find_normalization_4(5, 8, frames[i], fskymask)
        # temp_sky3 = find_normalization_4(9, 12, frames[i], fskymask)
        # temp_sky4 = find_normalization_4(13, 16, frames[i], fskymask)
        
        parent1, child1 = mp.Pipe(duplex=False)
        p1 = mp.Process(target=find_normalization_4, args=(1, 4, frames[i], fskymask, child1))
        p1.start()

        parent2, child2 = mp.Pipe(duplex=False)
        p2 = mp.Process(target=find_normalization_4, args=(5, 8, frames[i], fskymask, child2))
        p2.start()

        parent3, child3 = mp.Pipe(duplex=False)
        p3 = mp.Process(target=find_normalization_4, args=(9, 12, frames[i], fskymask, child3))
        p3.start()

        parent4, child4 = mp.Pipe(duplex=False)
        p4 = mp.Process(target=find_normalization_4, args=(13, 16, frames[i], fskymask, child4))
        p4.start()

        temp_sky1 = parent1.recv()
        temp_sky2 = parent2.recv()
        temp_sky3 = parent3.recv()
        temp_sky4 = parent4.recv()

        temp_sky = temp_sky1 + temp_sky2 + temp_sky3 + temp_sky4

        if i == n - 1:
            p1.join()
            p2.join()
            p3.join()
            p4.join()

        f.close()
        sg.close()

        # For each filter, there are 4 sky value, pick the largest one as the single normalization constant for that fhilter
        for k in range(4):
            sky_max[i][k] = max(temp_sky[k*4], temp_sky[k*4+1], temp_sky[k*4+2], temp_sky[k*4+3])
        print('Sky of '+frames[i]+':'+'\t'+str(sky_max[i]))

    hlist = []

    # hlist = [combine_nightflat(x, n, nx, ny, skymask[x].data, frames, sky_max) for x in range(1, 17)]
    parent1, child1 = mp.Pipe(duplex=False)
    p1 = mp.Process(target=combine_nightflat_4, args=(1, 4, n, nx, ny, fskymask, frames, sky_max, child1))
    p1.start()

    parent2, child2 = mp.Pipe(duplex=False)
    p2 = mp.Process(target=combine_nightflat_4, args=(5, 8, n, nx, ny, fskymask, frames, sky_max, child2))
    p2.start()

    parent3, child3 = mp.Pipe(duplex=False)
    p3 = mp.Process(target=combine_nightflat_4, args=(9, 12, n, nx, ny, fskymask, frames, sky_max, child3))
    p3.start()

    parent4, child4 = mp.Pipe(duplex=False)
    p4 = mp.Process(target=combine_nightflat_4, args=(13, 16, n, nx, ny, fskymask, frames, sky_max, child4))
    p4.start()

    hlist1 = parent1.recv()
    hlist2 = parent2.recv()
    hlist3 = parent3.recv()
    hlist4 = parent4.recv()

    p1.join()
    p2.join()
    p3.join()
    p4.join()

    hlist = hlist1 + hlist2 + hlist3 + hlist4

    hdu0 = fits.PrimaryHDU()
    hdu0.header = fits.getheader(filepath+frames[0])
    hlist.insert(0,hdu0)
    hduA = fits.HDUList(hlist)
    hduA.writeto(filepath+'night_flat_'+filt.lower()+'_'+month_str+'.fits', overwrite=True)
    skymask.close()


combine_frames(filt='ASU1')
combine_frames(filt='ASU2')
combine_frames(filt='ASU3')
combine_frames(filt='ASU4')




