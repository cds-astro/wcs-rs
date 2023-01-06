//! This module is an implementation of the SIP standard
//!
//! For the SIP convention, see
//! "The SIP convention for Representing Distortion in FITS Image Headers" by David L. Shupe et al.
//! in the proceedings of ADASS XIV (2005).

use fitsrs::hdu::Header;
use mapproj::sip::{Sip, SipAB, SipCoeff};

use crate::error::Error;

use crate::utils::string_to_keyword_type;

/// A method that return sip coefficients
fn retrieve_sip_coeffs(header: &Header, id: &'static str) -> Result<Option<SipCoeff>, Error> {
    let kw_order = format!("{}_ORDER ", id);
    let kw_order = unsafe { string_to_keyword_type(&kw_order) };
    if let Some(num_order) = header.get_parsed::<i64>(kw_order) {
        let num_order = num_order?;

        let coeffs = (0..=num_order)
            .flat_map(|i| (0..=num_order).map(move |j| (i, j)))
            .filter(|(i, j)| i + j <= num_order)
            .map(|(i, j)| {
                let kw_coeff_ij = format!("{}_{}_{}   ", id, i, j);
                let kw_coeff_ij = unsafe { string_to_keyword_type(&kw_coeff_ij) };

                header.get_parsed::<f64>(kw_coeff_ij).unwrap_or(Ok(0.0))
            })
            .collect::<Result<Vec<_>, fitsrs::error::Error>>()?
            .into_boxed_slice();

        Ok(Some(SipCoeff::new(coeffs)))
    } else {
        Ok(None)
    }
}

pub fn parse_sip(header: &Header, crpix1: f64, crpix2: f64) -> Result<Sip, Error> {
    // proj SIP coefficients
    let a_coeffs = retrieve_sip_coeffs(header, "A")?.unwrap_or_else(|| SipCoeff::new(Box::new([])));
    let b_coeffs = retrieve_sip_coeffs(header, "B")?.unwrap_or_else(|| SipCoeff::new(Box::new([])));

    let ap_coeffs = retrieve_sip_coeffs(header, "AP")?;
    let bp_coeffs = retrieve_sip_coeffs(header, "BP")?;

    let ab_proj = SipAB::new(a_coeffs, b_coeffs);

    let ab_deproj = match (ap_coeffs, bp_coeffs) {
        (Some(ap_coeffs), Some(bp_coeffs)) => Some(SipAB::new(ap_coeffs, bp_coeffs)),
        _ => None,
    };

    let naxis1 = *header.get_axis_size(1).ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS1"))? as f64;
    let naxis2 = *header.get_axis_size(2).ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS2"))? as f64;

    let u = (-crpix1)..=(naxis1 - crpix1);
    let v = (-crpix2)..=(naxis2 - crpix2);
    Ok(Sip::new(ab_proj, ab_deproj, u, v))
}
