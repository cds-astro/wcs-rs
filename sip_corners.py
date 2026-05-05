#!/usr/bin/env python3
"""
Compute the 4 corner (lon, lat) sky coordinates of a TAN-SIP FITS image.
"""

from astropy.wcs import WCS
from astropy.io import fits
import sys

def get_corners_lonlat(fits_path: str) -> list[tuple[float, float]]:
    with fits.open(fits_path) as hdul:
        # Find the first HDU with WCS/image data
        hdu = next(h for h in hdul if h.header.get("NAXIS1") and h.header.get("NAXIS2"))
        header = hdu.header

    wcs = WCS(header, naxis=2)

    naxis1 = header["NAXIS1"]  # width  (x, pixel columns)
    naxis2 = header["NAXIS2"]  # height (y, pixel rows)

    # FITS pixels are 1-indexed; corners in (x, y) = (col, row) order
    # (1,1) is bottom-left, (naxis1, naxis2) is top-right
    corners_pix = [
        (1,      1     ),  # bottom-left
        (naxis1, 1     ),  # bottom-right
        (naxis1, naxis2),  # top-right
        (1,      naxis2),  # top-left
    ]

    corners_sky = wcs.all_pix2world(corners_pix, 1)  # 1 = FITS 1-indexed origin

    return [(float(lon), float(lat)) for lon, lat in corners_sky]


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "image.fits"
    corners = get_corners_lonlat(path)
    labels = ["bottom-left", "bottom-right", "top-right", "top-left"]
    print(f"{'Corner':<15} {'Lon (deg)':>12} {'Lat (deg)':>12}")
    print("-" * 42)
    for label, (lon, lat) in zip(labels, corners):
        print(f"{label:<15} {lon:>12.6f} {lat:>12.6f}")
