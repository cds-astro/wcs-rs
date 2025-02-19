//! This module is an implementation of the SIP standard
//!
//! For the SIP convention, see
//! "The SIP convention for Representing Distortion in FITS Image Headers" by David L. Shupe et al.
//! in the proceedings of ADASS XIV (2005).

use mapproj::sip::{Sip, SipAB, SipCoeff};

use crate::error::Error;

use crate::params::WCSParams;

/// A method that return sip coefficients
macro_rules! combi_sip_coeff {
    ($params:ident, $id:tt, $( $p_q:literal ),+ ) => {
        paste! {
            Ok(vec![
                $(
                    $params.[< $id:lower _ $p_q >]
                ),*
            ])
        }
    };
}

use paste::paste;

macro_rules! build_sip_coeffs {
    ($params:ident, $id:tt) => {
        paste! {
            {
                let order = $params.[< $id:lower _ order >].unwrap_or(0);

                let coeffs: Vec<Option<f64>> = match order {
                    0 => combi_sip_coeff!(
                        $params,
                        $id,
                        0_0
                    ),
                    1 => combi_sip_coeff!(
                        $params,
                        $id,
                        0_0, 0_1,
                        1_0
                    ),
                    2 => combi_sip_coeff!(
                        $params,
                        $id,
                        0_0, 0_1, 0_2,
                        1_0, 1_1,
                        2_0
                    ),
                    3 => combi_sip_coeff!(
                        $params,
                        $id,
                        0_0, 0_1, 0_2, 0_3,
                        1_0, 1_1, 1_2,
                        2_0, 2_1,
                        3_0
                    ),
                    4 => combi_sip_coeff!(
                        $params,
                        $id,
                        0_0, 0_1, 0_2, 0_3, 0_4,
                        1_0, 1_1, 1_2, 1_3,
                        2_0, 2_1, 2_2,
                        3_0, 3_1,
                        4_0
                    ),
                    5 => combi_sip_coeff!(
                        $params,
                        $id,
                        0_0, 0_1, 0_2, 0_3, 0_4, 0_5,
                        1_0, 1_1, 1_2, 1_3, 1_4,
                        2_0, 2_1, 2_2, 2_3,
                        3_0, 3_1, 3_2,
                        4_0, 4_1,
                        5_0
                    ),
                    6 => combi_sip_coeff!(
                        $params,
                        $id,
                        0_0, 0_1, 0_2, 0_3, 0_4, 0_5, 0_6,
                        1_0, 1_1, 1_2, 1_3, 1_4, 1_5,
                        2_0, 2_1, 2_2, 2_3, 2_4,
                        3_0, 3_1, 3_2, 3_3,
                        4_0, 4_1, 4_2,
                        5_0, 5_1,
                        6_0
                    ),
                    _ => Err(crate::error::Error::SIPMaxOrderLimitReached)
                }?;

                Ok(
                    SipCoeff::new(
                        coeffs.into_iter()
                            .map(|c| c.unwrap_or(0.0))
                            .collect::<Vec<_>>()
                            .into_boxed_slice()
                    )
                )
            }
        }
    };
}

pub fn parse_sip(params: &WCSParams, naxis1: i64, naxis2: i64, crpix1: f64, crpix2: f64) -> Result<Sip, Error> {
    // proj SIP coefficients
    let a_coeffs: Result<SipCoeff, Error> = build_sip_coeffs!(params, A);
    let b_coeffs: Result<SipCoeff, Error> = build_sip_coeffs!(params, B);

    let ab_proj = SipAB::new(a_coeffs?, b_coeffs?);

    let ab_deproj = match (params.ap_order, params.bp_order) {
        (Some(_), Some(_)) => {
            let ap_coeffs: Result<SipCoeff, Error> = build_sip_coeffs!(params, AP);
            let bp_coeffs: Result<SipCoeff, Error> = build_sip_coeffs!(params, BP);

            Some(SipAB::new(ap_coeffs?, bp_coeffs?))
        }
        _ => None,
    };

    let naxis1 = naxis1 as f64;
    let naxis2 = naxis2 as f64;

    let u = (-crpix1)..=(naxis1 - crpix1);
    let v = (-crpix2)..=(naxis2 - crpix2);
    Ok(Sip::new(ab_proj, ab_deproj, u, v))
}
