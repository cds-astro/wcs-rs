//! This module deals with parsing the projection specific keywords
//!
//! Zenithal, pseudo-cylindrical, cylindrical and conic projections
//! are supported. Bonne, Polyconic and cubic projections are not

use std::f64::consts::PI;

use mapproj::{
    conic::{cod::Cod, coe::Coe, coo::Coo, cop::Cop},
    cylindrical::{car::Car, cea::Cea, cyp::Cyp, mer::Mer},
    hybrid::hpx::Hpx,
    math::HALF_PI,
    pseudocyl::{ait::Ait, mol::Mol, par::Par, sfl::Sfl},
    zenithal::{
        air::Air,
        arc::Arc,
        azp::Azp,
        ncp::Ncp,
        sin::{Sin, SinSlant},
        stg::Stg,
        szp::Szp,
        tan::Tan,
        zea::Zea,
        zpn::Zpn,
    },
    CanonicalProjection, CenteredProjection, LonLat,
};

use crate::params::WCSParams;

use super::Error;

#[derive(PartialEq, Debug)]
pub enum FiducialPoint {
    /// Fiducial point is the north pole (0, 90deg). This is the default for zenital projections
    NorthPole,
    /// Fiducial point is (0, 0). This is the default for cyl and pseudocyl projections
    Origin,
    UserSpeficied {
        /// Native Longitude of the fiducial point
        phi_0: f64,
        /// Native Latitude of the fiducial point
        theta_0: f64,
    },
}

/// Compute the celestial pole position
///
/// # Arguments
///
/// * `phi_p` - native longitude of the celestial pole. Corresponds to the LONPOLE fits keyword
/// * `theta_p` - native latitude of the celestial pole. Corresponds to the LATPOLE fits keyword
/// * `alpha_0` - celestial longitude of the fiducial point. Corresponds to the CRVAL1 fits keyword
/// * `delta_0` - celestial latitude of the fiducial point. Corresponds to the CRVAL2 fits keyword
/// * `fiducial_point` - the fiducial point position in the native system
///
/// # Returns
///
/// The celestial position of the pole
fn celestial_pole(
    phi_p: f64,
    theta_p: f64,
    alpha_0: f64,
    delta_0: f64,
    fiducial_point: &FiducialPoint,
) -> Result<LonLat, Error> {
    let (phi_0, theta_0) = match fiducial_point {
        FiducialPoint::NorthPole => {
            return Ok(LonLat::new(alpha_0, delta_0));
        }
        FiducialPoint::Origin => (0.0, 0.0),
        FiducialPoint::UserSpeficied { phi_0, theta_0 } => (*phi_0, *theta_0),
    };

    let (s_t0, c_t0) = theta_0.sin_cos();
    let (s_phi, c_phi) = (phi_p - phi_0).sin_cos();
    let (s_d0, c_d0) = delta_0.sin_cos();

    let a = c_t0 * c_t0 * s_phi * s_phi;
    let delta_p = if (a - 1.0).abs() < 1e-8 {
        if delta_0 == 0.0 {
            // Paper convention 6
            Ok(theta_p)
        } else {
            // Paper convention 3 error
            Err(Error::CelestialPoleInvalid)
        }
    } else {
        let b = (c_t0 * c_phi).atan2(s_t0);
        let c = (s_d0 / (1.0 - a).sqrt()).acos();

        let valid = -HALF_PI..=HALF_PI;

        match (valid.contains(&(b - c)), valid.contains(&(b + c))) {
            (false, false) => Err(Error::CelestialPoleInvalid),
            (true, false) => Ok(b - c),
            (false, true) => Ok(b + c),
            (true, true) => {
                // northerly solution chosen
                Ok(if theta_p >= 0.0 { b + c } else { b - c })
            }
        }
    }?;

    let alpha_p = if delta_0.abs() == HALF_PI {
        // Paper convention 1
        alpha_0
    } else if delta_p == HALF_PI {
        // Paper convention 2
        alpha_0 + phi_p - phi_0 - PI
    } else if delta_p == -HALF_PI {
        // Paper convention 2
        alpha_0 - phi_p + phi_0
    } else {
        let (s_dp, c_dp) = delta_p.sin_cos();

        let a = s_phi * c_t0 / c_d0;
        let b = (s_t0 - s_dp * s_d0) / (c_dp * c_d0);
        alpha_0 - a.atan2(b)
    };

    Ok(LonLat::new(alpha_p, delta_p))
}

pub trait WCSCanonicalProjection: CanonicalProjection {
    fn parse_proj(params: &WCSParams) -> Result<CenteredProjection<Self>, Error>
    where
        Self: Sized,
    {
        // Parse the celestial longitude of the fiducial point
        let crval1 = params.crval1.unwrap_or(0.0);
        // Parse the celestial latitude of the fiducial point
        let crval2 = params.crval2.unwrap_or(0.0);
        let crval: LonLat = LonLat::new(crval1.to_radians(), crval2.to_radians());

        let lonpole = params
            .lonpole
            .unwrap_or_else(|| if crval2 >= 0.0 { 0.0 } else { 180.0 });
        let latpole = params.latpole.unwrap_or(90.0);

        // Parse the native longitude of the fiducial point
        let native_fiducial_point = dbg!(match (params.pv1_1, params.pv1_2) {
            (Some(phi_0), Some(theta_0)) => FiducialPoint::UserSpeficied { phi_0, theta_0 },
            _ => Self::default_native_fiducial_point(params)?,
        });

        let positional_angle = if native_fiducial_point == FiducialPoint::NorthPole {
            -PI + lonpole.to_radians()
        } else {
            let pole = celestial_pole(
                lonpole.to_radians(),
                latpole.to_radians(),
                crval1.to_radians(),
                crval2.to_radians(),
                &native_fiducial_point,
            )?;

            // Compute the positional angle formed from (CRVAL1, CRVAL2) between (ALPHA_P, DELTA_P) and the north pole (0, 90°)
            let north_pole = LonLat::new(0.0, HALF_PI);

            let crval2pole_dist = crval.haversine_dist(&pole);
            let crval2np_dist = crval.haversine_dist(&north_pole);

            let pole2np_dist = pole.haversine_dist(&north_pole);

            let (s_02p, c_02p) = crval2pole_dist.sin_cos();
            let (s_02np, c_02np) = crval2np_dist.sin_cos();

            // A angle of a triangle on a sphere does not exceed PI
            // Use the law of cosines applied for geodesics on a sphere
            // https://www.theoremoftheday.org/GeometryAndTrigonometry/SphericalCos/TotDSphericalCos.pdf
            ((pole2np_dist.cos() - c_02p * c_02np) / (s_02p * s_02np)).acos()
        };

        let proj = Self::parse_internal_proj_params(params)?;

        let mut rotated_proj = CenteredProjection::new(proj);
        rotated_proj.set_proj_center_from_lonlat_and_positional_angle(&crval, positional_angle);
        dbg!(positional_angle.to_degrees());
        //rotated_proj.set_proj_center_from_lonlat(&crval);
        Ok(rotated_proj)
    }

    fn default_native_fiducial_point(params: &WCSParams) -> Result<FiducialPoint, Error>;

    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error>
    where
        Self: Sized;
}

// Zenithal projections
impl WCSCanonicalProjection for Azp {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // mu given in spherical radii, default value: 0.0
        let mu = params.pv2_1.unwrap_or(0.0);
        // gamma given in deg, default value: 0.0
        let gamma = params.pv2_2.unwrap_or(0.0);

        let azp = Azp::from_params(mu, gamma.to_radians());

        Ok(azp)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Szp {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // mu given in spherical radii, default value: 0.0
        let mu = params.pv2_1.unwrap_or(0.0);
        // phi_c given in deg, default value: 0.0
        let phi_c = params.pv2_2.unwrap_or(0.0);
        // theta_c given in deg, default value: 90.0
        let theta_c = params.pv2_3.unwrap_or(90.0);

        let szp = Szp::from_params(mu, phi_c.to_radians(), theta_c.to_radians());

        Ok(szp)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Tan {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Tan::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Stg {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Stg::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Sin {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Sin::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for SinSlant {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // xi dimensionless, default value: 0.0
        let xi = params.pv2_1.unwrap_or(0.0);
        // eta dimensionless, default value: 0.0
        let eta = params.pv2_2.unwrap_or(0.0);

        let sin_slant = SinSlant::new(-xi, eta);
        Ok(sin_slant)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Arc {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Arc::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Zpn {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        let mut coeffs = [
            params.pv2_0,
            params.pv2_1,
            params.pv2_2,
            params.pv2_3,
            params.pv2_4,
            params.pv2_5,
            params.pv2_6,
            params.pv2_7,
            params.pv2_8,
            params.pv2_9,
            params.pv2_10,
            params.pv2_11,
            params.pv2_12,
            params.pv2_13,
            params.pv2_14,
            params.pv2_15,
            params.pv2_16,
            params.pv2_17,
            params.pv2_18,
            params.pv2_19,
            params.pv2_20,
        ]
        .into_iter()
        .rev()
        .skip_while(|p| p.is_none())
        .map(|p| p.unwrap_or(0.0))
        .collect::<Vec<_>>();

        coeffs.reverse();

        Zpn::from_params(coeffs).ok_or(Error::InitProjection(
            Zpn::NAME,
            "negative polynomial in [0, pi]",
        ))
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Zea {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Zea::new())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Air {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // theta_b in deg, default value: 90.0
        let theta_b = params.pv2_1.unwrap_or(90.0);

        let airy = Air::from_param(theta_b.to_radians());
        Ok(airy)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

impl WCSCanonicalProjection for Ncp {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        let ncp = Ncp::new();
        Ok(ncp)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::NorthPole)
    }
}

// Cylindrical projections
impl WCSCanonicalProjection for Cyp {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // mu given in spherical radii, default value: 1.0
        let mu = params.pv2_1.unwrap_or(1.0);
        // lambda given in spherical radii, default value: 1.0
        let lambda = params.pv2_2.unwrap_or(1.0);

        let cyp = Cyp::from_params(mu, lambda);
        Ok(cyp)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Cea {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // lambda given in spherical radii, default value: 1.0
        let lambda = params.pv2_1.unwrap_or(1.0);

        let cea = Cea::from_param(lambda);
        Ok(cea)
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Car {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Car::default())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Mer {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Mer::default())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

// Pseudo-cylindrical projections
impl WCSCanonicalProjection for Sfl {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Sfl::default())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Par {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Par::default())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Mol {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Mol::default())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

impl WCSCanonicalProjection for Ait {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Ait::default())
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}

// Conic projections
impl WCSCanonicalProjection for Cop {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // theta_a given in deg, has no default value
        if let Some(theta_a) = params.pv2_1 {
            // eta given in deg, default value: 0
            let eta = params.pv2_2.unwrap_or(0.0);

            let cop = Cop::from_params(theta_a.to_radians(), eta.to_radians());
            Ok(cop)
        } else {
            Err(Error::InitProjection(
                Cop::NAME,
                "PV2_1 = theta_a must be defined as it has no default value",
            ))
        }
    }

    fn default_native_fiducial_point(params: &WCSParams) -> Result<FiducialPoint, Error> {
        if let Some(theta_a) = params.pv2_1 {
            Ok(FiducialPoint::UserSpeficied {
                phi_0: 0.0,
                theta_0: theta_a,
            })
        } else {
            Err(Error::InitProjection(
                Cop::NAME,
                "PV2_1 = theta_a must be defined as it has no default value",
            ))
        }
    }
}

impl WCSCanonicalProjection for Coe {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // theta_a given in deg, has no default value
        if let Some(theta_a) = params.pv2_1 {
            // eta given in deg, default value: 0
            let eta = params.pv2_2.unwrap_or(0.0);

            let cop = Coe::from_params(theta_a.to_radians(), eta.to_radians());
            Ok(cop)
        } else {
            Err(Error::InitProjection(
                Coe::NAME,
                "PV_1 = theta_a must be defined as it has no default value",
            ))
        }
    }

    fn default_native_fiducial_point(params: &WCSParams) -> Result<FiducialPoint, Error> {
        if let Some(theta_a) = params.pv2_1 {
            Ok(FiducialPoint::UserSpeficied {
                phi_0: 0.0,
                theta_0: theta_a,
            })
        } else {
            Err(Error::InitProjection(
                Coe::NAME,
                "PV2_1 = theta_a must be defined as it has no default value",
            ))
        }
    }
}

impl WCSCanonicalProjection for Cod {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // theta_a given in deg, has no default value
        if let Some(theta_a) = params.pv2_1 {
            // eta given in deg, default value: 0
            let eta = params.pv2_2.unwrap_or(0.0);

            let cod = Cod::from_params(theta_a.to_radians(), eta.to_radians());
            Ok(cod)
        } else {
            Err(Error::InitProjection(
                Cod::NAME,
                "PV_1 = theta_a must be defined as it has no default value",
            ))
        }
    }

    fn default_native_fiducial_point(params: &WCSParams) -> Result<FiducialPoint, Error> {
        if let Some(theta_a) = params.pv2_1 {
            Ok(FiducialPoint::UserSpeficied {
                phi_0: 0.0,
                theta_0: theta_a,
            })
        } else {
            Err(Error::InitProjection(
                Cod::NAME,
                "PV2_1 = theta_a must be defined as it has no default value",
            ))
        }
    }
}

impl WCSCanonicalProjection for Coo {
    fn parse_internal_proj_params(params: &WCSParams) -> Result<Self, Error> {
        // theta_a given in deg, has no default value
        if let Some(theta_a) = params.pv2_1 {
            // eta given in deg, default value: 0
            let eta = params.pv2_2.unwrap_or(0.0);

            let coo = Coo::from_params(theta_a.to_radians(), eta.to_radians());
            Ok(coo)
        } else {
            Err(Error::InitProjection(
                Coo::NAME,
                "PV_1 = theta_a must be defined as it has no default value",
            ))
        }
    }

    fn default_native_fiducial_point(params: &WCSParams) -> Result<FiducialPoint, Error> {
        if let Some(theta_a) = params.pv2_1 {
            Ok(FiducialPoint::UserSpeficied {
                phi_0: 0.0,
                theta_0: theta_a,
            })
        } else {
            Err(Error::InitProjection(
                Coo::NAME,
                "PV2_1 = theta_a must be defined as it has no default value",
            ))
        }
    }
}

impl WCSCanonicalProjection for Hpx {
    fn parse_internal_proj_params(_: &WCSParams) -> Result<Self, Error> {
        Ok(Hpx {})
    }

    fn default_native_fiducial_point(_: &WCSParams) -> Result<FiducialPoint, Error> {
        Ok(FiducialPoint::Origin)
    }
}
