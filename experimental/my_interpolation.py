from numpy import *
from scipy import interpolate
from astropy.io import fits



# LinearNDinterpolator

f1=f[1].data
c1=c[1].data

h=c1>0 # objects

f1[h] = nan

h2=~h # valid mask

coords = array(nonzero(h2)).T

values = f1[h2]

it = interpolate.LinearNDInterpolator(coords, values, fill_value=0)

nx = f1.shape[0]
ny = f1.shape[1]

x = linspace(0,nx-1, nx, dtype=int)
y = linspace(0,ny-1, ny, dtype=int)

x, y = meshgrid(x,y)

filled = it(x,y).T





# inter2d





