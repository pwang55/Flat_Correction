"""
Usage:

    In script folder:
    $ python weighted_combine_twilightflat.py path_to_files/flatlist.txt

    In corrected flats folder:
    $ python path_to_script/weighted_combine_twilightflat.py flatlist.txt

This script takes biassub_flat.corrected.fits and weighted average them to create twilight flats. You can have different filters/months in correctedlists. It will detect filter and month itself.

Input flatlist.txt should contain original raw flat file names.

"""
import numpy as np
from astropy.io import fits
import sys
import os


sig = 2.5


if len(sys.argv)!=2:
    print(__doc__)
    exit()

path_file = sys.argv[1]
filename = path_file.split('/')[-1]
if len(path_file.split('/'))>1:
    filepath = path_file[:-len(filename)]
else:
    filepath = ""

script_path = sys.path[0]
optic_dir = os.environ['_90PRIME_OPTIC_DIR']
bias_dark_flats_dir = os.environ['_90PRIME_BIAS_DARK_FLATS_DIR']


with open(path_file) as f:
    filelists = [ filepath+'biassub_'+l.strip().replace('.fits','.corrected.fits') for l in f ]

def weight_combine(filt='ASU1',month='Feb'):
    
    if month == 'Feb':
        mon = '02'
        month_dir = '/Feb/'
        month_str = 'feb'
    elif month == 'Mar':
        mon = '03'
        month_dir = '/Mar/'
        month_str = 'mar'

    # Grab flat files that have the correct filter and month
    flatfiles = [ l for l in filelists if (fits.getheader(l)['filter']==filt and fits.getheader(l)['date'].split('-')[1]==mon) ]

    # If there are no files in this filter/month, then don't execute the function
    if len(flatfiles) == 0:
        print('No files in '+filt+' and '+month+' , skip this')
        return

    n = len(flatfiles)

    fskymask = optic_dir+month_dir+'optic_flat_'+filt.lower()+'_skymask.fits'
    skymask = fits.open(fskymask)

    with fits.open(flatfiles[0]) as f:
        nx = f[1].header['naxis2']
        ny = f[1].header['naxis1']

    print('Using Skymask:\t'+month_dir+'optic_flat_'+filt.lower()+'_skymask.fits')
    print('\nCreating\t'+'flat_'+filt.lower()+'_'+month_str+'.fits using:')
    for x in range(n):
        print('\t'+flatfiles[x].split('/')[-1])
    print('\n')
    
    dat = np.zeros((n, nx, ny), dtype=float)
    dat2 = np.zeros((n, nx, ny), dtype=float)
    weights = np.zeros((n, nx, ny), dtype=float)

    hlist = []

    for j in range(1,17):
        for i in range(n):

            ho = skymask[j].data > 0.0

            print('Reading extension {} in {}'.format(j, flatfiles[i].split('/')[-1]))
            with fits.open(flatfiles[i]) as d:
                im = d[j].data
            dat[i] = im
            im[~ho] = 0.0
            weights[i][True] = np.median(im[im>0])/np.std(im[im>0])
            dat2[i] = dat[i]*weights[i]

        std_j = np.std(dat, axis=0)
        med_j = np.median(dat, axis=0)

        for i in range(n):
            hs = (dat[i] < med_j + sig*std_j) & (dat[i] > med_j - sig*std_j)
            weights[i][~hs] = 0.0
            dat2[i][~hs] = 0.0

        hduI = fits.ImageHDU()
        datf = np.sum(dat2, axis=0)/np.sum(weights, axis=0)
        h = np.isnan(datf)
        datf[h] = 0.0
        hduI.data = datf
        hduI.header = fits.getheader(flatfiles[0],j)
        hduI.header['BZERO'] = 0.
        hlist.append(hduI)

    hdu0 = fits.PrimaryHDU()
    hdu0.header = fits.getheader(flatfiles[0])
    hlist.insert(0,hdu0)
    hduA = fits.HDUList(hlist)
    hduA.writeto(filepath+'flat_{}_{}.fits'.format(filt.lower(),month_str))

    print('\nCreated: \t'+'flat_{}_{}.fits'.format(filt.lower(),month_str)+'\n')

    skymask.close()


weight_combine(filt='ASU1', month='Feb')
weight_combine(filt='ASU2', month='Feb')
weight_combine(filt='ASU3', month='Feb')
weight_combine(filt='ASU4', month='Feb')
weight_combine(filt='ASU1', month='Mar')
weight_combine(filt='ASU2', month='Mar')
weight_combine(filt='ASU3', month='Mar')
weight_combine(filt='ASU4', month='Mar')


