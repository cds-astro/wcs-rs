# WCS parsing library written in pure Rust

[![](https://img.shields.io/crates/v/wcs.svg)](https://crates.io/crates/wcs)
[![](https://img.shields.io/crates/d/wcs.svg)](https://crates.io/crates/wcs)
[![API Documentation on docs.rs](https://docs.rs/wcs/badge.svg)](https://docs.rs/wcs/)

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
* [X] Add hybrid projection (`HPX`)
* [X] Add pseudo cylindrical projections (`AIT`, `MOL`, `PAR`, `SFL`)
* [X] Add zenithal projections (`AIR`, `ARC`, `AZP`, `SIN`, `STG`, `SZP`, `TAN`, `ZEA`, `ZPN`)
* [ ] Zenithal projections not supported (`FEYE`, `NCP`)
* [ ] Add polyconic and pseudoconic projections (`BON, PCO`)?
* [ ] Add quad cube projections (`TSC`, `CSC`, `QSC`)?
* [ ] Check and possibly document constants to be added to match WCS projection bounds
* [X] Support `CRPIX` + `CD` convention
* [X] Support `CRPIX` + `PC` + `CDELT` convention
* [X] Support `CRPIX` + `CROTA` + `CDELT` convention
* [X] Add support for LONPOLE, (~LATPOLE approx)
* [X] SIP convention support
* [ ] SIP not tested
* [ ] [TPV](https://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html) World Coordinate System non-standard convention
* [ ] Add generation of projection files and plots (partially done)


Example
-------

```rust
use std::fs::File;
use std::io::BufReader;

use wcs::{
    WCS,
    WCSParams,
    ImgXY, LonLat
};

// Get the crval and crpix values along each axes
let crval1: f64 = 185.445488837;
let crval2: f64 = 4.47896032431;
let crpix1: f64 = 588.995094299;
let crpix2: f64 = 308.307905197;

// 2. Create a WCS from a specific header unit
let params: WCSParams = serde_json::from_str(r#"
    {
        "NAXIS": 2,
        "CTYPE1": "RA---TAN",
        "CTYPE2": "DEC--TAN",
        "EQUINOX": 2000.0,
        "LONPOLE": 180.0,
        "LATPOLE": 0.0,
        "CRVAL1": 185.445488837,
        "CRVAL2": 4.47896032431,
        "CRPIX1": 588.995094299,
        "CRPIX2": 308.307905197,
        "CUNIT1": "deg",
        "CUNIT2": "deg",
        "CD1_1": -0.000223666022989,
        "CD1_2": -0.000296578064584,
        "CD2_1": -0.000296427555509,
        "CD2_2": 0.000223774308964,
        "NAXIS1": 1080,
        "NAXIS2": 705
    }
"#).unwrap();
let wcs = WCS::new(&params).unwrap();

// 3. Once the WCS object is created, performs:
// * The projection of the center (lon, lat) = (crval1, crval2)
let lonlat = LonLat::new(crval1.to_radians(), crval2.to_radians());
let xy = wcs
    .proj_lonlat(&lonlat)
    .unwrap();
assert!((xy.x() - crpix1).abs() <= 1e-6);
assert!((xy.y() - crpix2).abs() <= 1e-6);

// * The unprojection of (X, Y) = (crpix1, crpix2)
let xy = ImgXY::new(crpix1, crpix2);
let lonlat = wcs
    .unproj_lonlat(&xy)
    .unwrap();
assert!((lonlat.lon() - crval1.to_radians()).abs() <= 1e-6);
assert!((lonlat.lat() - crval2.to_radians()).abs() <= 1e-6);

```