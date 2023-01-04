# WCS parsing library written in pure Rust

This crate must be used conjointly with [fitsrs](https://github.com/cds-astro/fitsrs)
It takes a FITS header and parse it to define a WCS.
A WCS object allow two operations:
* The projection of a (lon, lat) tuple given in degrees into the image space XY (in pixels)
* The unprojection of a vertex in the image space into a tuple (lon, lat) on the sphere
This crate heavily relies on the work done by F.-X. Pineau (@fxpineau) for the [mapproj](https://github.com/cds-astro/cds-mapproj-rust)
crate

The reference for the FITS WCS standard for undistorted celestial coordinate systems is Representations of celestial coordinates in FITS Calabretta, M. R., and Greisen, E. W., Astronomy & Astrophysics, 395, 1077-1122, 2002

To Do list
----------

* [X] Add conic projections (`COD`, `COE`, `COO`, `COP`)
* [X] Add cylindrical projections (`CAR`, `CEA`, `CYP`, `MER`)
* [ ] Add hybrid projection (`HPX`)
* [X] Add pseudo cylindrical projections (`AIT`, `MOL`, `PAR`, `SFL`)
* [X] Add zenithal projections (`AIR`, `ARC`, `AZP`, `SIN`, `STG`, `SZP`, `TAN`, `ZEA`, `ZPN`)
* [ ] Zenithal projections not supported (`FEYE`, `NCP`)
* [ ] Add polyconic and pseudoconic projections (`BON, PCO`)?
* [ ] Add quad cube projections (`TSC`, `CSC`, `QSC`)?
* [ ] Check and possibly document constants to be added to match WCS projection bounds
* [X] Support `CRPIX` + `CD` convention
* [X] Support `CRPIX` + `PC` + `CDELT` convention
* [X] Support `CRPIX` + `CROTA` + `CDELT` convention
* [ ] Add support for LONPOLE?
* [X] SIP convention support
* [ ] SIP not tested
* [ ] [TPV](https://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html) World Coordinate System non-standard convention
* [ ] Add generation of projection files and plots (partially done)


Example
-------

```rust
use std::fs::File;
use std::io::BufReader;

use fitsrs::Fits;
use wcsrs::{
    WCS,
    ImgXY, LonLat
};

fn main() {
    // 1. Parse a fits file using fitsrs
    let mut f = File::open("<your fits file path>").unwrap();

    let Fits { hdu: HDU { header, .. } } = Fits::from_reader(BufReader::new(f)).unwrap();
    // Get the crval and crpix values along each axes
    let crval1 = header
        .get_parsed::<f64>(b"CRVAL1  ")
        .unwrap_or(Ok(0.0))
        .unwrap();
    let crval2 = header
        .get_parsed::<f64>(b"CRVAL2  ")
        .unwrap_or(Ok(0.0))
        .unwrap();
    let crpix1 = header
        .get_parsed::<f64>(b"CRPIX1  ")
        .unwrap_or(Ok(0.0))
        .unwrap();
    let crpix2 = header
        .get_parsed::<f64>(b"CRPIX2  ")
        .unwrap_or(Ok(0.0))
        .unwrap();

    // 2. Create a WCS from a specific header unit
    let wcs = WCS::new(&header).unwrap();

    // 3. Once the WCS object is created, performs:
    // * The projection of the center (lon, lat) = (crval1, crval2)
    let lonlat = LonLat::new(crval1.to_radians(), crval2.to_radians());
    let xy = wcs
        .proj_lonlat(&lonlat)
        .unwrap();
    assert_delta!(xy.x(), crpix1 - 1.0, 1e-6);
    assert_delta!(xy.y(), crpix2 - 1.0, 1e-6);

    // * The unprojection of (X, Y) = (crpix1, crpix2)
    let xy = ImgXY::new(crpix1 - 1.0, crpix2 - 1.0);
    let lonlat = wcs
        .unproj_lonlat(&xy)
        .unwrap();
    assert_delta!(lonlat.lon(), crval1.to_radians(), 1e-6);
    assert_delta!(lonlat.lat(), crval2.to_radians(), 1e-6);
}
```