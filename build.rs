use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use paste::paste;
macro_rules! declare_card {
    (@ $keyword:ident, (($i:literal), ($( $j:literal ),*)), $type:ty) => {
        paste! {
            concat!(
                $(
                    "   ",
                    stringify!(
                        pub [< $keyword $i _ $j >]: $type,
                    ),
                    "\n"
                ),*
            )

        }
    };

    ($keyword:ident, (($( $i:literal ),*), $j:tt), $type:ty) => {
        concat!($(
            declare_card!(@ $keyword, (($i), $j), $type)
        ),*).as_bytes()
    };

    ($keyword:ident, ($( $i:literal ),*), $type:ty) => {
        paste! {
            concat!(
                $(
                    "   ",
                    stringify!(
                        pub [< $keyword $i >]: $type,
                    ),
                    "\n"
                ),*
            ).as_bytes()
        }
    };

    ($keyword:ident, $type:ty) => {
        concat!(
            "   ",
            stringify!(
                pub $keyword: $type,
            ),
            "\n"
        ).as_bytes()
    };
}

macro_rules! parse_card {
    (@@ ($header:ident, $key:ident, $optional:ident, $type:ty)) => {
        paste! {
            concat!(
                "           ",
                stringify!([<$key:lower>]),
                ": parse_",
                stringify!($optional),
                "_card_with_type::<",
                stringify!($type),
                ">(h, \"",
                stringify!($key),
                "\")?,",
                "\n"
            )
        }
    };

    (@ ($header:ident, $key:ident, $optional:ident, (), $type:ty)) => {
        parse_card!(@@ ($header, $key, $optional, $type))
    };

    (@ ($header:ident, $key:ident, $optional:ident, ($( $i:literal ),*), $type:ty)) => {
        paste! {
            concat!(
                $(
                    parse_card!(@@ ($header, [<$key $i>], $optional, $type))
                ),*
            )
        }
    };

    (@ ($header:ident, $key:ident, $optional:ident, [($( $i:literal ),*), $j:tt], $type:ty)) => {
        paste! {
            concat!(
                $(
                    parse_card!(@ ($header, [<$key $i _>], $optional, $j, $type))
                ),*
            )
        }
    };

    ($(($header:ident, $key:ident, $optional:ident, $i:tt, $type:ty)),*) => {
        concat!(
            $(
                parse_card!(@ ($header, $key, $optional, $i, $type))
            ),*
        ).as_bytes()
    };
}

pub fn write(path: PathBuf) -> Result<(), Box<dyn Error>> {
    let mut wcs_params_f = File::create(&path)?;
    wcs_params_f.write(
        r#"use fitsrs::{
    hdu::header::{extension::image::Image, Header},
};
use serde::{Deserialize, Serialize};
use fitsrs::card::CardValue;
use crate::{
    error::Error,
};

#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "UPPERCASE")]
pub struct WCSParams {
"#
        .as_bytes(),
    )?;
    // Mandatory keywords: NAXIS1, NAXIS2 and CTYPE1
    wcs_params_f.write(declare_card!(naxis, (1, 2), i64))?;
    wcs_params_f.write(declare_card!(ctype1, String))?;

    wcs_params_f.write(declare_card!(naxis, Option<i64>))?;
    wcs_params_f.write(declare_card!(naxis, (3, 4), Option<i64>))?;
    wcs_params_f.write(declare_card!(crpix, (1, 2, 3), Option<f64>))?;
    wcs_params_f.write(declare_card!(crval, (1, 2, 3), Option<f64>))?;
    wcs_params_f.write(declare_card!(cd, ((1, 2, 3), (1, 2, 3)), Option<f64>))?;
    wcs_params_f.write(declare_card!(pc, ((1, 2, 3), (1, 2, 3)), Option<f64>))?;
    wcs_params_f.write(declare_card!(cdelt, (1, 2, 3), Option<f64>))?;
    wcs_params_f.write(declare_card!(crota, (1, 2, 3), Option<f64>))?;
    wcs_params_f.write(declare_card!(ctype, (2, 3), Option<String>))?;
    wcs_params_f.write(declare_card!(pv1_, (0, 1, 2), Option<f64>))?;
    wcs_params_f.write(declare_card!(
        pv2_,
        (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20),
        Option<f64>
    ))?;
    wcs_params_f.write(declare_card!(epoch, Option<f64>))?;
    wcs_params_f.write(declare_card!(equinox, Option<f64>))?;
    wcs_params_f.write(declare_card!(radesys, Option<String>))?;
    wcs_params_f.write(declare_card!(lonpole, Option<f64>))?;
    wcs_params_f.write(declare_card!(latpole, Option<f64>))?;

    // SIP keywords
    // max order 6
    wcs_params_f.write(declare_card!(a_order, Option<i64>))?;

    wcs_params_f.write(declare_card!(
        a_,
        (0_0, 0_1, 0_2, 0_3, 0_4, 0_5, 0_6),
        Option<f64>
    ))?;
    wcs_params_f.write(declare_card!(
        a_,
        (1_0, 1_1, 1_2, 1_3, 1_4, 1_5),
        Option<f64>
    ))?;
    wcs_params_f.write(declare_card!(a_, (2_0, 2_1, 2_2, 2_3, 2_4), Option<f64>))?;
    wcs_params_f.write(declare_card!(a_, (3_0, 3_1, 3_2, 3_3), Option<f64>))?;
    wcs_params_f.write(declare_card!(a_, (4_0, 4_1, 4_2), Option<f64>))?;
    wcs_params_f.write(declare_card!(a_, (5_0, 5_1), Option<f64>))?;
    wcs_params_f.write(declare_card!(a_, (6_0), Option<f64>))?;

    wcs_params_f.write(declare_card!(b_order, Option<i64>))?;
    wcs_params_f.write(declare_card!(
        b_,
        (0_0, 0_1, 0_2, 0_3, 0_4, 0_5, 0_6),
        Option<f64>
    ))?;
    wcs_params_f.write(declare_card!(
        b_,
        (1_0, 1_1, 1_2, 1_3, 1_4, 1_5),
        Option<f64>
    ))?;
    wcs_params_f.write(declare_card!(b_, (2_0, 2_1, 2_2, 2_3, 2_4), Option<f64>))?;
    wcs_params_f.write(declare_card!(b_, (3_0, 3_1, 3_2, 3_3), Option<f64>))?;
    wcs_params_f.write(declare_card!(b_, (4_0, 4_1, 4_2), Option<f64>))?;
    wcs_params_f.write(declare_card!(b_, (5_0, 5_1), Option<f64>))?;
    wcs_params_f.write(declare_card!(b_, (6_0), Option<f64>))?;

    wcs_params_f.write(declare_card!(ap_order, Option<i64>))?;
    wcs_params_f.write(declare_card!(
        ap_,
        (0_0, 0_1, 0_2, 0_3, 0_4, 0_5, 0_6),
        Option<f64>
    ))?;
    wcs_params_f.write(declare_card!(
        ap_,
        (1_0, 1_1, 1_2, 1_3, 1_4, 1_5),
        Option<f64>
    ))?;
    wcs_params_f.write(declare_card!(ap_, (2_0, 2_1, 2_2, 2_3, 2_4), Option<f64>))?;
    wcs_params_f.write(declare_card!(ap_, (3_0, 3_1, 3_2, 3_3), Option<f64>))?;
    wcs_params_f.write(declare_card!(ap_, (4_0, 4_1, 4_2), Option<f64>))?;
    wcs_params_f.write(declare_card!(ap_, (5_0, 5_1), Option<f64>))?;
    wcs_params_f.write(declare_card!(ap_, (6_0), Option<f64>))?;

    wcs_params_f.write(declare_card!(bp_order, Option<i64>))?;
    wcs_params_f.write(declare_card!(
        bp_,
        (0_0, 0_1, 0_2, 0_3, 0_4, 0_5, 0_6),
        Option<f64>
    ))?;
    wcs_params_f.write(declare_card!(
        bp_,
        (1_0, 1_1, 1_2, 1_3, 1_4, 1_5),
        Option<f64>
    ))?;
    wcs_params_f.write(declare_card!(bp_, (2_0, 2_1, 2_2, 2_3, 2_4), Option<f64>))?;
    wcs_params_f.write(declare_card!(bp_, (3_0, 3_1, 3_2, 3_3), Option<f64>))?;
    wcs_params_f.write(declare_card!(bp_, (4_0, 4_1, 4_2), Option<f64>))?;
    wcs_params_f.write(declare_card!(bp_, (5_0, 5_1), Option<f64>))?;
    wcs_params_f.write(declare_card!(bp_, (6_0), Option<f64>))?;
    wcs_params_f.write(r#"
}

fn try_parse_card_from_header<T: CardValue>(header: &Header<Image>, key: &'static str) -> Result<Option<T>, Error> {
    let mut key_s = key.to_string();
    let filled_spaces = (0..(8 - key_s.len())).map(|_| " ").collect::<String>();
    key_s.push_str(&filled_spaces);
    
    let bytes = key_s.as_bytes().as_ptr() as *const [u8; 8];
    header.get_parsed::<T>(unsafe { &*bytes }).transpose().map_err(|e| e.into())
}
 
fn parse_optional_card_with_type<T: CardValue + std::str::FromStr>(header: &Header<Image>, key: &'static str) -> Result<Option<T>, Error> {
    match try_parse_card_from_header::<T>(header, key) {
        Ok(v) => Ok(v),
        _ => {
            let str = try_parse_card_from_header::<String>(header, key)
                .unwrap_or(None);

            Ok(if let Some(ss) = str {
                ss.trim().parse::<T>()
                    .map(|v| Some(v))
                    .unwrap_or(None)
                    //.map_err(|_| Error::CardWrongType(key.to_string(), std::any::type_name::<String>().to_string()))
            } else {
                // card not found but it is ok as it is not mandatory
                None
            })
        }
    }
}

fn parse_mandatory_card_with_type<T: CardValue>(header: &Header<Image>, key: &'static str) -> Result<T, Error> {
    match try_parse_card_from_header::<T>(header, key) {
        // No parsing error and found
        Ok(Some(v)) => {
            Ok(v)
        },
        // No error but not found, we return an error
        Ok(None) => Err(Error::MandatoryWCSKeywordsMissing(key)),
        // Return the parsing error
        Err(e) => {
            Err(e.into())
        }
    }
}

impl<'a> TryFrom<&'a Header<Image>> for WCSParams {
    type Error = Error;

    fn try_from(h: &'a Header<Image>) -> Result<Self, Self::Error> {
        Ok(WCSParams {"#.as_bytes())?;

    wcs_params_f.write(parse_card!(
        (h, NAXIS, mandatory, (1, 2), i64),
        (h, CTYPE1, mandatory, (), String),
        (h, CTYPE, optional, (2, 3), String),
        (h, NAXIS, optional, (), i64),
        (h, A_ORDER, optional, (), i64),
        (h, B_ORDER, optional, (), i64),
        (h, AP_ORDER, optional, (), i64),
        (h, BP_ORDER, optional, (), i64),
        (h, CRPIX, optional, (1, 2, 3), f64),
        (h, CRVAL, optional, (1, 2, 3), f64),
        (h, CROTA, optional, (1, 2, 3), f64),
        (h, CDELT, optional, (1, 2, 3), f64),
        (h, NAXIS, optional, (3, 4), i64),
        (h, LONPOLE, optional, (), f64),
        (h, LATPOLE, optional, (), f64),
        (h, EQUINOX, optional, (), f64),
        (h, EPOCH, optional, (), f64),
        (h, RADESYS, optional, (), String),
        (h, PV1_, optional, (0, 1, 2), f64),
        (
            h,
            PV2_,
            optional,
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20),
            f64
        ),
        (h, CD, optional, [(1, 2, 3), (1, 2, 3)], f64),
        (h, PC, optional, [(1, 2, 3), (1, 2, 3)], f64),
        (
            h,
            A_,
            optional,
            (
                0_0, 1_0, 2_0, 3_0, 4_0, 5_0, 6_0, 0_1, 1_1, 2_1, 3_1, 4_1, 5_1, 0_2, 1_2, 2_2,
                3_2, 4_2, 0_3, 1_3, 2_3, 3_3, 0_4, 1_4, 2_4, 0_5, 1_5, 0_6
            ),
            f64
        ),
        (
            h,
            AP_,
            optional,
            (
                0_0, 1_0, 2_0, 3_0, 4_0, 5_0, 6_0, 0_1, 1_1, 2_1, 3_1, 4_1, 5_1, 0_2, 1_2, 2_2,
                3_2, 4_2, 0_3, 1_3, 2_3, 3_3, 0_4, 1_4, 2_4, 0_5, 1_5, 0_6
            ),
            f64
        ),
        (
            h,
            B_,
            optional,
            (
                0_0, 1_0, 2_0, 3_0, 4_0, 5_0, 6_0, 0_1, 1_1, 2_1, 3_1, 4_1, 5_1, 0_2, 1_2, 2_2,
                3_2, 4_2, 0_3, 1_3, 2_3, 3_3, 0_4, 1_4, 2_4, 0_5, 1_5, 0_6
            ),
            f64
        ),
        (
            h,
            BP_,
            optional,
            (
                0_0, 1_0, 2_0, 3_0, 4_0, 5_0, 6_0, 0_1, 1_1, 2_1, 3_1, 4_1, 5_1, 0_2, 1_2, 2_2,
                3_2, 4_2, 0_3, 1_3, 2_3, 3_3, 0_4, 1_4, 2_4, 0_5, 1_5, 0_6
            ),
            f64
        )
    ))?;

    wcs_params_f.write(
        r#"
        })
    }
}
"#
        .as_bytes(),
    )?;

    Ok(())
}

fn main() {
    if let Err(err) = write("src/params.rs".into()) {
        // panic here for a nicer error message, otherwise it will
        // be flattened to one line for some reason
        panic!("Unable to generate shaders\n{}", err);
    }
}
