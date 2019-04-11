import cosmos
from pathlib import Path
from astropy.table import Table, hstack
import numpy as np
import astropy.units as u
import warnings
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.io.fits.verify import VerifyWarning


plt.close('all')


def load_reduced_data(filename):
    data = {}
    with fits.open(filename) as hdul:
        data['header'] = hdul[0].header
        data['sci'] = hdul['SCI'].data
        data['sources'] = Table.read(hdul['SEXTABLE'])
        data['wcs'] = WCS(hdul['SCI'].header)

    return data


def read_target_coordinates(header, *args, **kwargs):
    if not isinstance(header, fits.header.Header):
        header = fits.getheader(header, *args, **kwargs)

    ra = header['RA_TARG'] * u.deg
    dec = header['DEC_TARG'] * u.deg

    return SkyCoord(ra=ra, dec=dec, frame='icrs')


def get_targetindex(imageheader, catalog, nthneighbor=1):
    target_coordinates = read_target_coordinates(header)
    idx, sep2d, _ =  target_coordinates.match_to_catalog_sky(
                        cosmos.cat_to_sc(cosmoscat), nthneighbor=nthneighbor)

    return np.asscalar(idx), np.asscalar(sep2d)


def update_fitsfile(filename, extension, data, header=None):

    hdumaker = fits.BinTableHDU
    if not isinstance(data, Table):
        hdumaker = fits.ImageHDU

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', VerifyWarning)
        hdu = hdumaker(data, header=header, name=extension)

    with fits.open(filename, mode='update') as hdul:

        while True:
            try:
                hdul.__delitem__(extension)
            except KeyError:
                break

        hdul.append(hdu)
        hdul.flush()

# %matplotlib tk

# TODO: matching distance threshold related to psf size
# TODO: matching maybe with respect to magnitude. Search around 2d and then
# selection condition

# Set up directory and filenames

# no sources: 3, 6
mapidx = 10
catdir = Path().home() / 'data/catalogs'
mapdir = Path().home() / 'data/reduced_maps/dthresh_3'

cosmoscatname = 'UVISTA_final_v4.1.fits'
sourcefilename = 'file{:02d}_backsubTrue.fits'.format(mapidx)  # 'thirdmap_withskysubtraction_drz.fits'

# Load data
cosmoscat = Table.read(catdir / cosmoscatname)
cosmoscat.meta['NAME'] = 'UVISTA4.1'
data = load_reduced_data(mapdir / sourcefilename)
wcs = data['wcs']
sources = data['sources']
image = data['sci']
header = data['header']
pixscale = proj_plane_pixel_scales(wcs) * u.deg
targetidx, _ = get_targetindex(header, cosmoscat)

MATCHRADIUS = 2 * u.arcsec  # 0.2 * u.arcsec
RACORRECTION = 0.187 * u.arcsec  # - 0.05 * u.arcsec
DECCORRECTION = .3 * u.arcsec  # 0.2 * u.arcsec
# %matplotlib tk

sources.keys()


# autocorrection function that takes median of dec and ra deviation

# ax = cosmos.scatter_catalog_in_image(cosmoscat, wcs, image, scatterall=True)
# cosmos.scatter_catalog_in_image(sources, wcs)

sources['Y_WORLD'] = sources['Y_WORLD'] + DECCORRECTION.to_value(u.deg)
sources['X_WORLD'] = sources['X_WORLD'] + RACORRECTION.to_value(u.deg)

sources['X_WORLD_MATCHCORRECTED'] = sources['X_WORLD'].copy()
sources['Y_WORLD_MATCHCORRECTED'] = sources['Y_WORLD'].copy()

# % % match sources
with warnings.catch_warnings():
    warnings.simplefilter('ignore', UserWarning)
    sources, sep2d = cosmos.match_sources(sources, cosmoscat, dist_threshold=None)
    sources, sep2d = cosmos.match_sources(sources, cosmoscat, dist_threshold=MATCHRADIUS)

# sources.masked
# % %
sources['Y_WORLD'] = sources['Y_WORLD'] - DECCORRECTION.to_value(u.deg)
sources['X_WORLD'] = sources['X_WORLD'] - RACORRECTION.to_value(u.deg)

sources['sep2d'] = sep2d

h = fits.header.Header()
h['REFCAT'] = 'uvista'
h['matchrad'] = MATCHRADIUS.to_value(u.deg), 'matchradius in deg'
h['RACORR'] = RACORRECTION.to_value(u.deg), 'RA correction in degrees for match'
h['DECCORR'] = DECCORRECTION.to_value(u.deg), 'DEC correction in degrees for match'

# update_fitsfile(mapdir / sourcefilename, 'MATCHEDCAT', sources, h)
# % %
#%matplotlib tk

cosmos.mask_double_matches(sources, cosmoscat)
matchidx = sources[cosmoscat.meta['NAME']]
nmatches = np.sum(~matchidx.mask)
print('found {:d} matches for {} sources ({:.2f})'.format(
                            nmatches, len(sources), nmatches/len(sources)))
#  % %

sources_with_match, matched_counterpart = cosmos.delete_unmatched_rows(
                                                            sources, cosmoscat)

cosmos.plot_match_coordinate_deviation(sources_with_match, matched_counterpart)

# sources_with_match['Y_WORLD'] = sources_with_match['Y_WORLD'] - DECCORRECTION.to_value(u.deg)
# sources_with_match['X_WORLD'] = sources_with_match['X_WORLD'] - RACORRECTION.to_value(u.deg)
# %matplotlib tk


# % % plot detected and matched sources

assert len(sources_with_match) == len(matched_counterpart)

ax = cosmos.scatter_catalog_in_image(sources_with_match, wcs=wcs,
                                     image=image, scatterall=True, color='g',
                                     label='detected and cosmos counterpart')
cosmos.scatter_catalog_in_image(matched_counterpart, wcs=wcs,
                                scatterall=True, ax=ax, color='r',
                                label='cosmos match')
cosmos.connect_catalog_points(sources_with_match, matched_counterpart, wcs=wcs,
                              ax=ax)

ax.legend()

# %%


# for colname in matched_counterpart.colnames:
#     matched_counterpart[colname].name = matched_counterpart[colname].name + '_refcat'
# matchedtable = hstack([sources_with_match, matched_counterpart])




# header catalog [uvista, mccracken], matchradius

# %%
'''
# %% Plot Restframe Color
rf_UmV = Table.read(catdir / 'UVISTA_final_v4.1_153-155rf.fits')
rf_VmJ = Table.read(catdir / 'UVISTA_final_v4.1_155-161rf.fits')

assert len(rf_UmV) == len(rf_VmJ)
assert len(cosmoscat) == len(rf_UmV)

VmJ = -2.5 * np.log10(rf_VmJ['L155']/rf_VmJ['L161'])
UmV = -2.5 * np.log10(rf_UmV['L153']/rf_UmV['L155'])

# % % Plot mass

M05_cat = Table.read(catdir / 'UVISTA_final_M05_v4.1_fout.fits')
mass = M05_cat['lmass']
redshift = M05_cat['z']
redshift[targetidx]
len(M05_cat)
'''
