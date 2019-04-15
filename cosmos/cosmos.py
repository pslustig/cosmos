from astropy.table import Table, MaskedColumn
from astropy.coordinates import SkyCoord, match_coordinates_sky
import warnings
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from .zscale import zscale


def get_coordinate_column_names(cat):

    namepaires = (('ra', 'dec'), ('RA', 'DEC'), ('X_WORLD', 'Y_WORLD'))

    for raname, decname in namepaires:
        if (raname in cat.keys()) and (decname in cat.keys()):
            return raname, decname

    raise ValueError('Coordinate columns not found. Available keys '
                     'are\n{}'.format(cat.keys()))


def delete_unmatched_rows(catalog, reference, return_referenceorder=False):
    catalog = catalog[~catalog[reference.meta['NAME']].mask]
    reference = reference[catalog[reference.meta['NAME']]]

    ret = catalog, reference
    if return_referenceorder:
        ret = catalog, reference, catalog[reference.meta['NAME']]

    return ret


def find_double_matches(catalog, matchcolname):
    matchidx = catalog[matchcolname]
    unique, inv, counts = np.unique(matchidx, return_inverse=True,
                                    return_counts=True)
    # get idx in unique array of sources with multiple matches
    multipleidx = np.arange(len(counts))[counts > 1]
    # get idx of lines in catalog that have multiple matches
    doublelines = np.arange(len(catalog))[np.isin(inv, multipleidx)]

    return doublelines, unique[counts > 1], counts[counts > 1]


def mask_double_matches(catalog, reference=None, matchcolname=None,
                        returnmaskonly=False):

    if matchcolname is None:
        matchcolname = reference.meta['NAME']
    doublelines, _, _ = find_double_matches(catalog, matchcolname)

    if not returnmaskonly:
        catalog[matchcolname].mask[doublelines] = True
    return doublelines


def print_double_matches(catalog, reference):
    doublelines, doublevalue, counts = find_double_matches(catalog, reference)
    print('sources {} ({} in catalog) matched more than one time ({})'.format(
          doublevalue.tolist(), doublelines, counts))


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

    cols = get_coordinate_column_names(cat)
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


def plot_match_coordinate_deviation(catalog1, catalog2):
    '''
    catalog1 and 2 must provide ra and dec columns with ra and dec in degrees
    '''
    # TODO: make coordinates first
    try:
        coords1 = SkyCoord(catalog1['X_WORLD_MATCHCORRECTED'],
                           catalog1['Y_WORLD_MATCHCORRECTED'], unit='deg')
        print('Using corrected coordinates for matchplot')
    except KeyError:
        coords1 = cat_to_sc(catalog1)

    coords2 = cat_to_sc(catalog2)
    x = (coords1.ra - coords2.ra).to_value(u.arcsec)
    x *= np.cos(coords1.dec.to_value(u.rad))
    y = (coords1.dec - coords2.dec).to_value(u.arcsec)

    fig, ax = plt.subplots(1, 2, sharex=True, sharey=True)
    ax[0].scatter(x, y, alpha=.3)
    ax[1].hist2d(x, y, bins=100)

    fs = 15
    ax[0].set_xlabel('(ra1-ra1) * cos(dec1) [arcsec]', fontsize=fs)
    ax[1].set_xlabel('(ra1-ra1) * cos(dec1) [arcsec]', fontsize=fs)
    ax[0].set_ylabel('dec1-dec2 [arcsec]', fontsize=fs)

    return ax


def match_sources(catalog, reference, dist_threshold):

    catalog_sc = cat_to_sc(catalog)
    reference_sc = cat_to_sc(reference)
    idx, sep2d, _ = match_coordinates_sky(catalog_sc, reference_sc)

    # catalog = Table(catalog, masked=True)
    mask = np.zeros(len(sep2d), dtype=bool)
    if dist_threshold is not None:
        mask = sep2d > dist_threshold

    catalog[reference.meta['NAME']] = MaskedColumn(idx, mask=mask)
    return catalog, sep2d


def catalog_world_to_pix(catalog, wcs):

    cols = get_coordinate_column_names(catalog)
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
        vmin, vmax = zscale(image)
        ax.imshow(image, origin='lower', cmap='Greys', vmin=vmin, vmax=vmax)

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
