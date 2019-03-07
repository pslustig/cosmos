from astropy.table import Table
from pathlib import Path
import astropy.units as u
from astropy.io.fits.header import Header
import warnings
from astropy.io import fits

warnings.warn('EVERYTHING IS HARDCODED FOR UVISTA_final_M05_v4.1.fout FILE')

catdir = Path().home() / 'data/catalogs'
catalog = Table.read(catdir / 'UVISTA_final_M05_v4.1.fout', format='ascii')
print('catalog loaded')
# %%

colnames = ('id', 'z', 'ltau', 'metal', 'lage', 'Av', 'lmass', 'lsfr',
            'lssfr', 'la2t', 'chi2')
colunits = (None, None, u.dex(u.yr), None, u.dex(u.yr), None, u.dex(u.M_sun),
            u.dex(u.M_sun/u.yr), u.dex(u.yr), None, None)

colunits = (None, None, u.dex(u.yr), None, u.dex(u.yr), None, u.dex(u.M_sun),
            u.dex(u.M_sun/u.yr), u.dex(u.yr), 'log(age/tau)', None)


for column, newname in zip(catalog.colnames, colnames):
    catalog.rename_column(column, newname)

header = Header()
header['FAST_VER'] = '0.9b', 'Fast Version'
header['PHOT_CAT'] = ('UVISTA_final_v4.1.cat', 'Photometric catalog file')
header['Z_CAT'] = 'UVISTA_final_v4.1.zout', 'Photometric redshift file'
header['TEMP_ERR'] = 'TEMPLATE_ERROR.fast.v0.2', 'Template error function'
header['AB_ZP'] = 25.00
header['LIBRARY'] = 'Maraston (2005)'
header['SFH'] = 'Exponentially declining SFH: SFR ~ exp(-t/tau)'
header['IMF'] = 'Kroupa', 'Stellar IMF'
header['METALIC'] = 0.020, 'metallicity'

header['LTAUYMIN'] = 7, 'minimum log(tau/yr)'
header['LTAUYMAX'] = 10, 'maximum log(tau/yr)'
header['LTAUYSTP'] = 0.10, 'log(tau/yr) steps'

header['LAGEYMIN'] = 7, 'minimum log(age/yr)'
header['LAGEYMAX'] = 10.1, 'maximum log(age/yr)'
header['LAGEYSTP'] = 0.10, 'log(age/yr) steps'

header['A_VMIN'] = 0.0, 'minimum A_V'
header['A_VMAX'] = 4.0, 'maximum A_V'
header['A_VSTP'] = 0.10, 'A_V steps'

header['ZMIN'] = 0.01, 'minimum z'
header['ZMAX'] = 6.00, 'maximum z'
header['ZSTP'] = 0.01, 'z steps'


for colname, colunit in zip(colnames, colunits):
    if colunit is not None:
        s = colunit
        if not isinstance(colunit, str):
            s = colunit.to_string()
        header['{}UNIT'.format(colname)] = s

filters = ('226', '225', '224', '223', '21', '20', '19', '18', '83', '82',
           '81', '79', '80', '78', '88', '184', '186', '190', '192', '194',
           '195', '181', '183', '185', '188', '193', '197', '120', '121')

for i, filter in enumerate(filters):
    header['FILTER{}'.format(i)] = int(filter)

phdu = fits.PrimaryHDU(header=header)
thdu = fits.BinTableHDU(data=catalog)
# %%
fits.HDUList([phdu, thdu]).writeto(catdir / 'UVISTA_final_M05_v4.1_fout.fits')
