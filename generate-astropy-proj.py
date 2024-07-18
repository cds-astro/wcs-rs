from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import glob
import numpy as np

input_files = glob.glob('./examples/*.fits')
print(input_files)

for filename in input_files:
    print(filename)
    # Load the FITS hdulist using astropy.io.fits
    hdulist = fits.open(filename)

    # Parse the WCS keywords in the primary HDU
    w = WCS(hdulist[0].header)

    naxis1 = hdulist[0].header.get("NAXIS1")
    naxis2 = hdulist[0].header.get("NAXIS2")
    naxis3 = hdulist[0].header.get("NAXIS3")

    # Get random sky coordinates
    X = np.random.random(100) * naxis1
    Y = np.random.random(100) * naxis2
    Z = np.random.random(100) * 0

    try:
        coord = w.wcs_pix2world(np.vstack((X, Y)).T, 1)
        ra, dec = ((coord[:, 0] * u.deg).to_value(u.rad), (coord[:, 1] * u.deg).to_value(u.rad))
    except:
        coord = w.wcs_pix2world(np.vstack((X, Y, Z)).T, 1)
        ra, dec = ((coord[:, 0] * u.deg).to_value(u.rad), (coord[:, 1] * u.deg).to_value(u.rad))

    tab = np.column_stack((ra, dec, X, Y))
    np.savetxt(filename + ".csv", tab, delimiter=",")
