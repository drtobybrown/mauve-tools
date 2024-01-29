# Foreground dust extinction for the MAUVE / VERTICO galaxies
from astropy.table import Table
from astropy.coordinates import SkyCoord
from dustmaps.config import config
import dustmaps.sfd
import numpy as np

# ensure the dust maps are downloaded
config['data_dir'] = '/arc/projects/mauve/bin/dustmaps/'
dustmaps.sfd.fetch()

# read in VERTICO galaxies
df = Table.read('/arc/projects/vertico/share/catalogues/VERTICO_basic.fits')
ebmv = np.zeros(len(df))
coords = SkyCoord(ra = df["RA"], dec=df["DEC"], unit='deg', frame='fk5')

# look up the EBmV for each position
sfd = dustmaps.sfd.SFDQuery()
for gal, coord in zip(df["Galaxy"], coords):
    ebmv = sfd(coord)
    print("{gal} EBmV = {ebmv}".format(gal=gal, ebmv=ebmv))
