import cosmos
from pathlib import Path
from astropy.table import Table
import numpy as np
import astropy.units as u
import warnings
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales


# Set up directory and filenames

# %matplotlib tk
catdir = Path().home() / 'data/catalogs'
mapdir = Path().home() / 'data/reduced_maps/dthresh_3'

cosmoscatname = 'UVISTA_final_v4.1.fits'
sourcefilename = 'file02_backsubFalse.fits'  # 'thirdmap_withskysubtraction_drz.fits'

rf_UmV = Table.read(catdir / 'UVISTA_final_v4.1_153-155rf.fits')
rf_VmJ = Table.read(catdir / 'UVISTA_final_v4.1_155-161rf.fits')

assert len(rf_UmV) == len(rf_VmJ)



rf_VmJ.keys()
rf_UmV.keys()




x = -2.5 * np.log10(rf_VmJ['L155']/rf_VmJ['L161'])
y = -2.5 * np.log10(rf_UmV['L153']/rf_UmV['L155'])

select = rf_VmJ['z'] < .1
x = x[select]
y = y[select]

fig, ax = plt.subplots(1, 1)
ax.scatter(x, y, alpha=.2)
