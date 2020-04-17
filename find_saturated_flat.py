"""

Usage:

    In script directory:
    $ python find_saturated_flat.py path_to_file/flatlist.txt

    In files directory:
    $ python path_to_script/find_saturated_flat.py flatlist.txt

Default saturated level is 65500 for the 99 percentile. This script will return names of flat files that are saturated.


"""
from numpy import *
from astropy.io import fits
from glob import glob
import sys


path_file = sys.argv[1]
filename = path_file.split('/')[-1]
if len(path_file.split('/'))>1:
    filepath = path_file[:-len(filename)]
else:
    filepath = ""

filelists = []
with open(path_file) as f:
    filelists = [ l.strip() for l in f ]


for x in range(len(filelists)):

    f = fits.open(filepath+filelists[x])

    for j in range(1,17):
            dat = f[j].data
	    top = percentile(dat,99)
	    if top > 65500.0:
		print filelists[x]
		break

    f.close()



