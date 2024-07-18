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

impl WCSParams {
    fn from(
        sky_center: (f64, f64),
        image_size: (f64, f64),
        angular_res: (f64, f64),
        projection: &str,
        coo_system: CooSystem,
        positional_angle: f64,
    ) -> Self {
        todo!();
    }
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
            naxis: h.get_parsed::<i64>(b"NAXIS   ").transpose()?,
            naxis1,
            naxis2,
            crpix1: h.get_parsed::<f64>(b"CRPIX1  ").transpose()?,
            crpix2: h.get_parsed::<f64>(b"CRPIX2  ").transpose()?,
            crval1: h.get_parsed::<f64>(b"CRVAL1  ").transpose()?,
            crval2: h.get_parsed::<f64>(b"CRVAL2  ").transpose()?,
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
            equinox: h
                .get_parsed::<f64>(b"EQUINOX ")
                .transpose()
                .unwrap_or_else(|_| {
                    let s = h
                        .get_parsed::<String>(b"EQUINOX ")
                        .transpose()
                        .unwrap_or_else(|_| None);

                    if let Some(s) = s {
                        s.parse::<f64>().ok()
                    } else {
                        None
                    }
                }),
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
