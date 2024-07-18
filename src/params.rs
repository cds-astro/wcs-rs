use fitsrs::{
    card::CardValue,
    hdu::header::{extension::image::Image, Header},
};
use serde::{Deserialize, Serialize};

use crate::{
    coo_system::{self, CooSystem},
    error::Error,
    utils,
};

use paste::paste;

#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "UPPERCASE")]
pub struct WCSParams {
    /// number of axis
    pub naxis: Option<i64>,
    /// size of the first axis
    pub naxis1: u64,
    /// size of the second axis
    pub naxis2: u64,
    pub crpix1: Option<f64>,
    pub crpix2: Option<f64>,
    pub crval1: Option<f64>,
    pub crval2: Option<f64>,
    /// rotation in the pixel space
    pub cd1_1: Option<f64>,
    pub cd1_2: Option<f64>,
    pub cd2_1: Option<f64>,
    pub cd2_2: Option<f64>,

    pub pc1_1: Option<f64>,
    pub pc1_2: Option<f64>,
    pub pc2_1: Option<f64>,
    pub pc2_2: Option<f64>,

    pub cdelt1: Option<f64>,
    pub cdelt2: Option<f64>,

    pub crota2: Option<f64>,

    pub ctype1: String,
    pub ctype2: Option<String>,

    pub epoch: Option<f64>, /// deprecated in favor of EQUINOX
    pub equinox: Option<f64>,
    pub radesys: Option<String>,

    pub pv1_0: Option<f64>,
    /// native longitude fiducial point
    pub pv1_1: Option<f64>,
    /// native latitude fiducial point
    pub pv1_2: Option<f64>,

    /// native longitude of the celestial pole
    pub lonpole: Option<f64>,
    /// native latitude of the celestial pole
    pub latpole: Option<f64>,

    /// projections parameters
    pub pv2_0: Option<f64>,
    pub pv2_1: Option<f64>,
    pub pv2_2: Option<f64>,
    pub pv2_3: Option<f64>,
    pub pv2_4: Option<f64>,
    pub pv2_5: Option<f64>,
    pub pv2_6: Option<f64>,
    pub pv2_7: Option<f64>,
    pub pv2_8: Option<f64>,
    pub pv2_9: Option<f64>,
    pub pv2_10: Option<f64>,
    pub pv2_11: Option<f64>,
    pub pv2_12: Option<f64>,
    pub pv2_13: Option<f64>,
    pub pv2_14: Option<f64>,
    pub pv2_15: Option<f64>,
    pub pv2_16: Option<f64>,
    pub pv2_17: Option<f64>,
    pub pv2_18: Option<f64>,
    pub pv2_19: Option<f64>,
    pub pv2_20: Option<f64>,
}


use std::f64::consts::PI;

impl WCSParams {
    /// lon and lat in degrees, equivalent to CRVAL1 and CRVAL2
    /// lon and lat in degrees of the celestial pole
    /// in pixels
    /// in deg/px, equivalent to CDELT1 and CDELT2
    /// projection
    /// coo system e.g. ICRS, GALACTIC, ..
    fn from<P: mapproj::Projection>(
        celestial_ref: (f64, f64),
        celestial_pole: (f64, f64),
        pixel_size: (f64, f64),
        angular_res: (f64, f64),
        proj: &P,
        coo_system: CooSystem,
    ) -> Result<Self, Error> {
        /*let proj_name = proj.short_name();

        let lonpole = match proj_name {
            "AIT" | "MOL" | "MER" | "CAR" => {
                
            }
            // zenithal
            "SIN" | "TAN" => {
                (PI - positional_angle.to_radians()).to_degrees()
            }
            // conic
            "COD" => {

            }
            // hybrid
            "HPX" => {}
            _ => return Error::NotImplementedProjection(proj_name.to_string())
        }*/
        todo!();
    }
}


macro_rules! generate_wcs_param {
    ($param:ident, $( $i:literal )+, $t:ty) => {
        paste! {
            pub struct WCSa {
                $(
                    [<$param $i>]: $t
                ),*
            }
        }
    }
}

generate_wcs_param!(A, 10, f64);


macro_rules! parse_non_mandatory_card_with_type {
    ($header:ident, $e:expr, $t:ty) => {
        {
            let key = match $e.len() {
                1 => concat!($e, "       "),
                2 => concat!($e, "      "),
                3 => concat!($e, "     "),
                4 => concat!($e, "    "),
                5 => concat!($e, "   "),
                6 => concat!($e, "  "),
                7 => concat!($e, " "),
                8 => $e,
                _ => unreachable!()
            };
            let bytes = key.as_bytes().as_ptr() as *const [u8; 8];

            $header.get_parsed::<$t>(unsafe { &*bytes }).transpose()
        }
    };
    ($header:ident, $e:expr, $t:ty, $( $ts:ty ),+) => {
        match parse_non_mandatory_card_with_type!($header, $e, $t) {
            Ok(v) => Ok(v),
            _ => {
                let value = parse_non_mandatory_card_with_type!($header, $e, $( $ts )*)?;

                if let Some(value) = value {
                    value.parse::<$t>()
                        .map(|v| Some(v))
                        .map_err(|_| Error::CardWrongType($e.to_string(), std::any::type_name::<$t>().to_string()))
                } else {
                    // card not found but it is ok as it is not mandatory
                    Ok(None)
                }
            }
        }
    };
}

impl<'a> TryFrom<&'a Header<Image>> for WCSParams {
    type Error = Error;

    fn try_from(h: &'a Header<Image>) -> Result<Self, Self::Error> {
        let xtension = h.get_xtension();

        let naxis1 = *xtension
            .get_naxisn(1)
            .ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS1"))?;
        let naxis2 = *xtension
            .get_naxisn(2)
            .ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS2"))?;

        Ok(WCSParams {
            naxis: parse_non_mandatory_card_with_type!(h, "NAXIS", i64)?,
            naxis1,
            naxis2,
            crpix1: parse_non_mandatory_card_with_type!(h, "CRPIX1", f64)?,
            crpix2: parse_non_mandatory_card_with_type!(h, "CRPIX2", f64)?,
            crval1: parse_non_mandatory_card_with_type!(h, "CRVAL1", f64)?,
            crval2: parse_non_mandatory_card_with_type!(h, "CRVAL2", f64)?,
            cd1_1: h.get_parsed::<f64>(b"CD1_1   ").transpose()?,
            cd1_2: h.get_parsed::<f64>(b"CD1_2   ").transpose()?,
            cd2_1: h.get_parsed::<f64>(b"CD2_1   ").transpose()?,
            cd2_2: h.get_parsed::<f64>(b"CD2_2   ").transpose()?,
            pc1_1: h.get_parsed::<f64>(b"PC1_1   ").transpose()?,
            pc1_2: h.get_parsed::<f64>(b"PC1_2   ").transpose()?,
            pc2_1: h.get_parsed::<f64>(b"PC2_1   ").transpose()?,
            pc2_2: h.get_parsed::<f64>(b"PC2_2   ").transpose()?,
            cdelt1: h.get_parsed::<f64>(b"CDELT1  ").transpose()?,
            cdelt2: h.get_parsed::<f64>(b"CDELT2  ").transpose()?,
            crota2: h.get_parsed::<f64>(b"CROTA2  ").transpose()?,
            ctype1: utils::retrieve_mandatory_parsed_keyword(h, "CTYPE1  ")?,
            ctype2: h.get_parsed::<String>(b"CTYPE2  ").transpose()?,
            epoch: h.get_parsed::<f64>(b"EPOCH   ").transpose()?,
            equinox: parse_non_mandatory_card_with_type!(h, "EQUINOX", f64, String)?,
            radesys: h.get_parsed::<String>(b"RADESYS ").transpose()?,
            pv1_0: h.get_parsed::<f64>(b"PV1_0   ").transpose()?,
            pv1_1: h.get_parsed::<f64>(b"PV1_1   ").transpose()?,
            pv1_2: h.get_parsed::<f64>(b"PV1_2   ").transpose()?,
            lonpole: h.get_parsed::<f64>(b"LONPOLE ").transpose()?,
            latpole: h.get_parsed::<f64>(b"LATPOLE ").transpose()?,
            pv2_0: h.get_parsed::<f64>(b"PV2_0   ").transpose()?,
            pv2_1: h.get_parsed::<f64>(b"PV2_1   ").transpose()?,
            pv2_2: h.get_parsed::<f64>(b"PV2_2   ").transpose()?,
            pv2_3: h.get_parsed::<f64>(b"PV2_3   ").transpose()?,
            pv2_4: h.get_parsed::<f64>(b"PV2_4   ").transpose()?,
            pv2_5: h.get_parsed::<f64>(b"PV2_5   ").transpose()?,
            pv2_6: h.get_parsed::<f64>(b"PV2_6   ").transpose()?,
            pv2_7: h.get_parsed::<f64>(b"PV2_7   ").transpose()?,
            pv2_8: h.get_parsed::<f64>(b"PV2_8   ").transpose()?,
            pv2_9: h.get_parsed::<f64>(b"PV2_9   ").transpose()?,
            pv2_10: h.get_parsed::<f64>(b"PV2_10  ").transpose()?,
            pv2_11: h.get_parsed::<f64>(b"PV2_11  ").transpose()?,
            pv2_12: h.get_parsed::<f64>(b"PV2_12  ").transpose()?,
            pv2_13: h.get_parsed::<f64>(b"PV2_13  ").transpose()?,
            pv2_14: h.get_parsed::<f64>(b"PV2_14  ").transpose()?,
            pv2_15: h.get_parsed::<f64>(b"PV2_15  ").transpose()?,
            pv2_16: h.get_parsed::<f64>(b"PV2_16  ").transpose()?,
            pv2_17: h.get_parsed::<f64>(b"PV2_17  ").transpose()?,
            pv2_18: h.get_parsed::<f64>(b"PV2_18  ").transpose()?,
            pv2_19: h.get_parsed::<f64>(b"PV2_19  ").transpose()?,
            pv2_20: h.get_parsed::<f64>(b"PV2_20  ").transpose()?,
        })
    }
}
