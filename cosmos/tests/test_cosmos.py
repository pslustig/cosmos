from __future__ import absolute_import, division, print_function

# from .. import cosmos relative import failes here, do not know why
import pytest
import cosmos
import numpy as np
from astropy.table import Table
import astropy.units as u


def test_pos_in_shape():
    pos = np.arange(16).reshape((8, 2)) - 4
    shape = (5, 5)
    inside = cosmos.position_in_shape(pos, shape)
    assert np.array_equal(pos[2:4], inside)


def test_pos_in_shape_extra():
    pos = np.arange(16).reshape((8, 2)) - 4
    extra = np.arange(len(pos))
    shape = (5, 5)
    inside, inside_extra = cosmos.position_in_shape(pos, shape, extra=extra)
    assert np.array_equal(pos[2:4], inside)
    assert np.array_equal(inside_extra, np.arange(2, 4))


def test_cat_to_sc():
    cat1 = Table()
    cat1['ra'] = [1, 5] * u.deg
    cat1['dec'] = [1, 5] * u.deg

    cosmos.cat_to_sc(cat1)
    cat2 = Table()
    cat2['X_WORLD'] = [1.1, 4.9] * u.deg
    cat2['Y_WORLD'] = [.9, 5.1] * u.deg
    cosmos.cat_to_sc(cat2)


def test_cat_to_sc_warnings():
    cat1 = Table()
    cat1['ra'] = [1, 5]
    cat1['dec'] = [1, 5]

    with pytest.warns(UserWarning):
        cosmos.cat_to_sc(cat1)


def test_match_sources():
    cat1 = Table()
    cat1['ra'] = [1, 5] * u.deg
    cat1['dec'] = [1, 5] * u.deg
    cat2 = Table()
    cat2['X_WORLD'] = [1.1, 4.9] * u.deg
    cat2['Y_WORLD'] = [.9, 5.1] * u.deg
    cat2.meta['NAME'] = 'to_match'

    cosmos.match_sources(cat1, cat2, dist_threshold=1*u.deg)
    assert np.array_equal(cat1['to_match'], np.arange(2))


def test_match_sources_mask():
    cat1 = Table()
    cat1['ra'] = [1, 5] * u.deg
    cat1['dec'] = [1, 5] * u.deg
    cat2 = Table()
    cat2['X_WORLD'] = [1.1, 4.9] * u.deg
    cat2['Y_WORLD'] = [.9, 5.1] * u.deg
    cat2.meta['NAME'] = 'to_match'
    cosmos.match_sources(cat1, cat2, dist_threshold=.1*u.deg)
    assert np.array_equal(cat1['to_match'].mask, [True, True])
