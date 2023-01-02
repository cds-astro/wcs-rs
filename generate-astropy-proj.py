from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import glob
import numpy as np

input_files = glob.glob('./examples/*.fits')
print(input_files)

for filename in input_files:
    # Load the FITS hdulist using astropy.io.fits
    hdulist = fits.open(filename)

    # Parse the WCS keywords in the primary HDU
    w = WCS(hdulist[0].header)

    naxis1 = hdulist[0].header.get("NAXIS1")
    naxis2 = hdulist[0].header.get("NAXIS2")
    
    # Get random sky coordinates
    X = np.random.random(100) * naxis1
    Y = np.random.random(100) * naxis2

    coord = w.pixel_to_world(X, Y)
    ra, dec = (coord.ra.to_value(u.rad), coord.dec.to_value(u.rad))

    tab = np.column_stack((ra, dec, X, Y))
    np.savetxt(filename + ".csv", tab, delimiter=",")
