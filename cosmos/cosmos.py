from astropy.table import MaskedColumn
from astropy.coordinates import SkyCoord, match_coordinates_sky
import warnings
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt


def cat_to_sc(cat):
    """Extract positions from cat and return corresponding SkyCoord
    Parameters
    ----------
    cat : :class:`astropy.table.Table`
        a table containing sky columns with units
    Returns
    -------
    :class:`astropy.coordinates.SkyCoord`
        the corresponding SkyCoord object
    """

    if 'ra' in cat.keys() and 'dec' in cat.keys():
        cols = ['ra', 'dec']
    elif 'X_WORLD' in cat.keys() and 'Y_WORLD' in cat.keys():
        cols = ['X_WORLD', 'Y_WORLD']

    ra_unit, dec_unit = cat[cols[0]].unit, cat[cols[1]].unit
    if ra_unit is None:
        warnings.warn('ra unit is not provided, assuming degrees')
        ra_unit = u.deg
    if dec_unit is None:
        warnings.warn('dec unit is not provided, assuming degrees')
        dec_unit = u.deg

    coords = SkyCoord(cat[cols[0]], cat[cols[1]],
                      unit=(ra_unit, dec_unit))

    return coords


def match_sources(catalog, reference, dist_threshold):

    catalog_sc = cat_to_sc(catalog)
    reference_sc = cat_to_sc(reference)
    idx, sep2d, _ = match_coordinates_sky(catalog_sc, reference_sc)
    mask = sep2d > dist_threshold
    catalog[reference.meta['NAME']] = MaskedColumn(idx, mask=mask)
    return sep2d


def catalog_world_to_pix(catalog, wcs):

    if 'ra' in catalog.keys() and 'dec' in catalog.keys():
        cols = ['ra', 'dec']
    elif 'X_WORLD' in catalog.keys() and 'Y_WORLD' in catalog.keys():
        cols = ['X_WORLD', 'Y_WORLD']

    pos = np.array(wcs.wcs_world2pix(catalog[cols[0]], catalog[cols[1]], 0)).T

    return pos


def position_in_shape(pos, shape, within=(0, 1), extra=None):
    limits = shape * np.asarray(within)[:, np.newaxis]
    inside = np.sum((pos >= limits[0]) & (pos <= limits[1] - 1), 1) == 2

    pos = pos[inside]
    if extra is not None:
        pos = (pos, extra[inside])

    return pos


def scatter_catalog_in_image(catalog, wcs, image=None, ax=None,
                             scatterall=False, **kwargs):

    shape = wcs._naxis
    pospix = catalog_world_to_pix(catalog, wcs)
    if not scatterall:
        pospix = position_in_shape(
                            pospix, shape=shape)

    if ax is None:
        fig, ax = plt.subplots(1, 1)

    if image is not None:
        ax.imshow(image, origin='lower', cmap='Greys')

    ax.scatter(pospix[:, 0], pospix[:, 1], **kwargs)

    ax.set_xlim(0, shape[0])
    ax.set_ylim(0, shape[1])

    return ax


def connect_catalog_points(catalog1, catalog2, wcs, ax=None):
    assert len(catalog1) == len(catalog2), 'catalogs must have the same length'
    cat1pos = catalog_world_to_pix(catalog1, wcs)
    cat2pos = catalog_world_to_pix(catalog2, wcs)

    if ax is None:
        fig, ax = plt.subplots(1, 1)

    for c1, c2 in zip(cat1pos, cat2pos):
        ax.plot([c1[0], c2[0]], [c1[1], c2[1]])
