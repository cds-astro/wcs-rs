extern crate mapproj;
#[macro_use]
extern crate quick_error;

pub mod error;

use coo_system::CooSystem;
use error::Error;
pub mod coo_system;
pub mod params;
mod projection;
mod sip;
mod utils;

use crate::projection::WCSCanonicalProjection;
pub use params::WCSParams;

// Imports
use mapproj::{
    conic::{cod::Cod, coe::Coe, coo::Coo, cop::Cop},
    cylindrical::{car::Car, cea::Cea, cyp::Cyp, mer::Mer},
    hybrid::hpx::Hpx,
    img2celestial::Img2Celestial,
    img2proj::{WcsImgXY2ProjXY, WcsWithSipImgXY2ProjXY},
    pseudocyl::{ait::Ait, mol::Mol, par::Par, sfl::Sfl},
    zenithal::{
        air::Air, arc::Arc, azp::Azp, ncp::Ncp, sin::Sin, stg::Stg, szp::Szp, tan::Tan, zea::Zea,
        zpn::Zpn,
    },
    XYZ,
};

use paste::paste;
/// macro
macro_rules! create_specific_proj {
    ( $proj_name:ident, $params:expr, $ctype1:expr, $naxis1:expr, $naxis2:expr, $crpix1:expr, $crpix2:expr, $img2proj:expr ) => {{
        let (positional_angle, proj) = $proj_name::parse_proj(&$params)?;

        let is_sip_found = &$ctype1[($ctype1.len() - 3)..] == "SIP";
        if is_sip_found {
            let sip = sip::parse_sip($params, $naxis1, $naxis2, $crpix1, $crpix2)?;
            let img2proj = WcsWithSipImgXY2ProjXY::new($img2proj, sip);

            paste! {
                Ok((WCSCelestialProj::[ <$proj_name Sip> ](Img2Celestial::new(img2proj, proj)), positional_angle))
            }
        } else {
            Ok((
                WCSCelestialProj::$proj_name(Img2Celestial::new($img2proj, proj)),
                positional_angle,
            ))
        }
    }};
}

/// Structure alias coming from mapproj defining
/// image space pixel coordinates
pub type ImgXY = mapproj::ImgXY;
/// Structure alias coming from mapproj defining
/// longitude and latitude expressed in degrees
pub type LonLat = mapproj::LonLat;

#[derive(Debug)]
pub struct WCS {
    /* Metadata keywords */
    /// Size of the image in pixels in its i-th dimension
    naxisi: Box<[i64]>,
    /// Field of view of the image along NAXIS1
    fov1: f64,
    /// Field of view of the image along NAXIS2
    fov2: f64,
    /// Main sub structure defining the projection
    proj: WCSProj,
}

/// Main object structure descripting a WCS object
/// Once created, the user can proceed two operation on it
/// * The projection of a (lon, lat) tuple onto the image space.
///   Results are given in pixels
/// * The unprojection of a (x, y) tuple given in pixel coordinates onto the sphere.
///   Results are given as a (lon, lat) tuple expressed in degrees
impl WCS {
    pub fn new(params: &WCSParams) -> Result<Self, Error> {
        let naxisi = match params.naxis {
            2 => {
                let naxis1 = params.naxis1.ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS1"))?;
                let naxis2 = params.naxis2.ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS2"))?;

                Ok(vec![naxis1, naxis2].into_boxed_slice())
            }
            3 => {
                let naxis1 = params.naxis1.ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS1"))?;
                let naxis2 = params.naxis2.ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS2"))?;
                let naxis3 = params.naxis3.ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS3"))?;

                Ok(vec![naxis1, naxis2, naxis3].into_boxed_slice())
            }
            4 => {
                let naxis1 = params.naxis1.ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS1"))?;
                let naxis2 = params.naxis2.ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS2"))?;
                let naxis3 = params.naxis3.ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS3"))?;
                let naxis4 = params.naxis4.ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS4"))?;

                Ok(vec![naxis1, naxis2, naxis3, naxis4].into_boxed_slice())
            }
            _ => {
                Err(Error::NotSupportedNaxis(params.naxis))
            }
        }?;

        // At least NAXIS >= 2
        let proj = WCSProj::new(naxisi[0], naxisi[1], params)?;

        let fov1 = proj.s_lon * (naxisi[0] as f64);
        let fov2 = proj.s_lat * (naxisi[1] as f64);

        Ok(WCS {
            naxisi,
            fov1,
            fov2,
            proj,
        })
    }

    /// Returns the dimensions of the image given by the NAXIS1 x NAXIS2 keyword
    pub fn img_dimensions(&self) -> &[i64] {
        &self.naxisi[..]
    }

    pub fn field_of_view(&self) -> (f64, f64) {
        (self.fov1, self.fov2)
    }

    /// Project a (lon, lat) in ICRS
    ///
    /// The result is given a (X, Y) tuple expressed in pixel coordinates.
    ///
    /// # Arguments
    ///
    /// * `lonlat`: the 3D sphere vertex expressed as a (lon, lat) tuple given in degrees
    pub fn proj(&self, lonlat: &LonLat) -> Option<ImgXY> {
        self.proj.proj_lonlat(lonlat)
    }

    /// Unproject a (X, Y) point to get a position on the sky in ICRS system
    ///
    /// # Arguments
    ///
    /// * `img_pos`: the image space point expressed as a (X, Y) tuple given en pixels
    pub fn unproj(&self, img_pos: &ImgXY) -> Option<LonLat> {
        self.proj.unproj_lonlat(img_pos)
    }

    /// Get the coordinate system frame
    pub fn coo_system(&self) -> &CooSystem {
        self.proj.coo_system()
    }
}

use std::ops::Deref;
impl Deref for WCS {
    type Target = WCSProj;

    fn deref(&self) -> &Self::Target {
        &self.proj
    }
}

pub struct WCSProj {
    /// The right part of the CTYPE keyword
    /// The projection type
    proj: WCSCelestialProj,
    /// The left part of the CTYPE keyword
    /// The coordinate system
    coo_system: CooSystem,
    pos_angle: f64,
    s_lon: f64, // scale in degrees along the longitude axis
    s_lat: f64, // scale in degrees along the latitude axis
}

impl std::fmt::Debug for WCSProj {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("WCS")
            .field("coosys", &self.coo_system)
            .field("euler angle", &self.pos_angle)
            .finish()
    }
}

/// Main enum structure descripting a WCS object
/// Once created, the user can proceed two operation on it
/// * The projection of a (lon, lat) tuple onto the image space.
///   Results are given in pixels
/// * The unprojection of a (x, y) tuple given in pixel coordinates onto the sphere.
///   Results are given as a (lon, lat) tuple expressed in degrees
pub enum WCSCelestialProj {
    // Zenithal
    Azp(Img2Celestial<Azp, WcsImgXY2ProjXY>),
    Szp(Img2Celestial<Szp, WcsImgXY2ProjXY>),
    Tan(Img2Celestial<Tan, WcsImgXY2ProjXY>),
    Stg(Img2Celestial<Stg, WcsImgXY2ProjXY>),
    Sin(Img2Celestial<Sin, WcsImgXY2ProjXY>),
    Arc(Img2Celestial<Arc, WcsImgXY2ProjXY>),
    Zpn(Img2Celestial<Zpn, WcsImgXY2ProjXY>),
    Zea(Img2Celestial<Zea, WcsImgXY2ProjXY>),
    Air(Img2Celestial<Air, WcsImgXY2ProjXY>),
    Ncp(Img2Celestial<Ncp, WcsImgXY2ProjXY>),
    // Cylindrical
    Cyp(Img2Celestial<Cyp, WcsImgXY2ProjXY>),
    Cea(Img2Celestial<Cea, WcsImgXY2ProjXY>),
    Car(Img2Celestial<Car, WcsImgXY2ProjXY>),
    Mer(Img2Celestial<Mer, WcsImgXY2ProjXY>),
    // Pseudo-Cylindrical
    Sfl(Img2Celestial<Sfl, WcsImgXY2ProjXY>),
    Par(Img2Celestial<Par, WcsImgXY2ProjXY>),
    Mol(Img2Celestial<Mol, WcsImgXY2ProjXY>),
    Ait(Img2Celestial<Ait, WcsImgXY2ProjXY>),
    // Conic
    Cop(Img2Celestial<Cop, WcsImgXY2ProjXY>),
    Cod(Img2Celestial<Cod, WcsImgXY2ProjXY>),
    Coe(Img2Celestial<Coe, WcsImgXY2ProjXY>),
    Coo(Img2Celestial<Coo, WcsImgXY2ProjXY>),
    // Hybrid
    Hpx(Img2Celestial<Hpx, WcsImgXY2ProjXY>),

    // SIP handling
    // Zenithal
    AzpSip(Img2Celestial<Azp, WcsWithSipImgXY2ProjXY>),
    SzpSip(Img2Celestial<Szp, WcsWithSipImgXY2ProjXY>),
    TanSip(Img2Celestial<Tan, WcsWithSipImgXY2ProjXY>),
    StgSip(Img2Celestial<Stg, WcsWithSipImgXY2ProjXY>),
    SinSip(Img2Celestial<Sin, WcsWithSipImgXY2ProjXY>),
    ArcSip(Img2Celestial<Arc, WcsWithSipImgXY2ProjXY>),
    ZpnSip(Img2Celestial<Zpn, WcsWithSipImgXY2ProjXY>),
    ZeaSip(Img2Celestial<Zea, WcsWithSipImgXY2ProjXY>),
    AirSip(Img2Celestial<Air, WcsWithSipImgXY2ProjXY>),
    NcpSip(Img2Celestial<Ncp, WcsWithSipImgXY2ProjXY>),
    // Cylindrical
    CypSip(Img2Celestial<Cyp, WcsWithSipImgXY2ProjXY>),
    CeaSip(Img2Celestial<Cea, WcsWithSipImgXY2ProjXY>),
    CarSip(Img2Celestial<Car, WcsWithSipImgXY2ProjXY>),
    MerSip(Img2Celestial<Mer, WcsWithSipImgXY2ProjXY>),
    // Pseudo-Cylindrical
    SflSip(Img2Celestial<Sfl, WcsWithSipImgXY2ProjXY>),
    ParSip(Img2Celestial<Par, WcsWithSipImgXY2ProjXY>),
    MolSip(Img2Celestial<Mol, WcsWithSipImgXY2ProjXY>),
    AitSip(Img2Celestial<Ait, WcsWithSipImgXY2ProjXY>),
    // Conic
    CopSip(Img2Celestial<Cop, WcsWithSipImgXY2ProjXY>),
    CodSip(Img2Celestial<Cod, WcsWithSipImgXY2ProjXY>),
    CoeSip(Img2Celestial<Coe, WcsWithSipImgXY2ProjXY>),
    CooSip(Img2Celestial<Coo, WcsWithSipImgXY2ProjXY>),
    // Hybrid
    HpxSip(Img2Celestial<Hpx, WcsWithSipImgXY2ProjXY>),
}

fn parse_pc_matrix(params: &WCSParams) -> Option<(f64, f64, f64, f64)> {
    let pc11 = params.pc1_1;
    let pc12 = params.pc1_2;
    let pc21 = params.pc2_1;
    let pc22 = params.pc2_2;

    let pc_matrix_found = match (&pc11, &pc12, &pc21, &pc22) {
        (None, None, None, None) => false,
        // The CD1_1 keyword has been found
        // We are in a case where the CDij are given
        _ => true,
    };

    if pc_matrix_found {
        let pc11 = pc11.unwrap_or(1.0);
        let pc12 = pc12.unwrap_or(0.0);
        let pc21 = pc21.unwrap_or(0.0);
        let pc22 = pc22.unwrap_or(1.0);

        Some((pc11, pc12, pc21, pc22))
    } else {
        None
    }
}

fn parse_cd_matrix(params: &WCSParams) -> Option<(f64, f64, f64, f64)> {
    let cd11 = params.cd1_1;
    let cd12 = params.cd1_2;
    let cd21 = params.cd2_1;
    let cd22 = params.cd2_2;

    let cd_matrix_found = match (&cd11, &cd12, &cd21, &cd22) {
        (None, None, None, None) => false,
        // The CD1_1 keyword has been found
        // We are in a case where the CDij are given
        _ => true,
    };

    if cd_matrix_found {
        let cd11 = cd11.unwrap_or(1.0);
        let cd12 = cd12.unwrap_or(0.0);
        let cd21 = cd21.unwrap_or(0.0);
        let cd22 = cd22.unwrap_or(1.0);

        Some((cd11, cd12, cd21, cd22))
    } else {
        None
    }
}

impl WCSProj {
    /// Create a WCS from a specific fits header parsed with fitsrs
    /// # Param
    /// * `naxis1` - Size of the image in its first dimension (in pixels)
    /// * `naxis2` - Size of the image in its second dimension (in pixels)
    /// * `params` - Header unit coming from fitsrs.
    ///   This contains all the cards of one HDU.
    pub fn new(naxis1: i64, naxis2: i64, params: &WCSParams) -> Result<Self, Error> {
        // 1. Identify the image <-> intermediate projection
        // a. Linear transformation matrix cases:
        // - CRPIXi + CDij
        // - CRPIXi + CDELTi + CROTA2
        // - CRPIXi + CDELTi + PCij
        let crpix1 = params.crpix1.unwrap_or(0.0);
        let crpix2 = params.crpix2.unwrap_or(0.0);

        // Choice of the wcs order:
        // 1 - Priority to define the projection is given to CD
        // 2 - Then, to the couple PC + CDELT
        // 3 - Finally to the old CROTA + CDELT convention
        let (img2proj, s_lon, s_lat) =
            if let Some((cd11, cd12, cd21, cd22)) = parse_cd_matrix(params) {
                // CDij case
                (
                    WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22),
                    cd11.abs(),
                    cd22.abs(),
                )
            } else {
                // Search for CDELTi
                let cdelt1 = params.cdelt1.unwrap_or(1.0);
                let cdelt2 = params.cdelt2.unwrap_or(1.0);

                if let Some((pc11, pc12, pc21, pc22)) = parse_pc_matrix(params) {
                    // CDELTi + PCij case
                    (
                        WcsImgXY2ProjXY::from_pc(
                            crpix1, crpix2, pc11, pc12, pc21, pc22, cdelt1, cdelt2,
                        ),
                        (cdelt1 * pc11).abs(),
                        (cdelt2 * pc22).abs(),
                    )
                } else {
                    // CDELTi + CROTA2 case
                    let crota2 = params.crota2.unwrap_or(0.0);
                    let cosc = crota2.to_radians().cos();

                    (
                        WcsImgXY2ProjXY::from_cr(crpix1, crpix2, crota2, cdelt1, cdelt2),
                        (cdelt1 * cosc).abs(),
                        (cdelt2 * cosc).abs(),
                    )
                }
            };

        // 2. Identify the projection type
        let ctype1 = &params.ctype1;
        let proj_name = &ctype1[5..=7];

        let (proj, pos_angle) = match proj_name.as_bytes() {
            // Zenithal
            b"AZP" => {
                create_specific_proj!(Azp, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"SZP" => {
                create_specific_proj!(Szp, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"TAN" => {
                create_specific_proj!(Tan, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"STG" => {
                create_specific_proj!(Stg, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"SIN" => {
                create_specific_proj!(Sin, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"ARC" => {
                create_specific_proj!(Arc, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"ZPN" => {
                create_specific_proj!(Zpn, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"ZEA" => {
                create_specific_proj!(Zea, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"AIR" => {
                create_specific_proj!(Air, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"NCP" => {
                create_specific_proj!(Ncp, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            // Cylindrical
            b"CYP" => {
                create_specific_proj!(Cyp, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"CEA" => {
                create_specific_proj!(Cea, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"CAR" => {
                create_specific_proj!(Car, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"MER" => {
                create_specific_proj!(Mer, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            // Pseudo-cylindrical
            b"SFL" => {
                create_specific_proj!(Sfl, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"PAR" => {
                create_specific_proj!(Par, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"MOL" => {
                create_specific_proj!(Mol, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"AIT" => {
                create_specific_proj!(Ait, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            // Conic
            b"COP" => {
                create_specific_proj!(Cop, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"COD" => {
                create_specific_proj!(Cod, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"COE" => {
                create_specific_proj!(Coe, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            b"COO" => {
                create_specific_proj!(Coo, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            // HEALPix
            b"HPX" => {
                create_specific_proj!(Hpx, params, ctype1, naxis1, naxis2, crpix1, crpix2, img2proj)
            }
            _ => Err(Error::NotImplementedProjection(proj_name.to_string())),
        }?;

        let coo_system = CooSystem::parse(&params)?;

        Ok(WCSProj {
            proj,
            coo_system,
            pos_angle,
            s_lon,
            s_lat,
        })
    }

    /// Project a (lon, lat) given in ICRS frame to get its corresponding location on the image
    ///
    /// The result is given a (X, Y) tuple expressed in pixel coordinates.
    ///
    /// # Arguments
    ///
    /// * `lonlat`: a coo expressed as (lon, lat) tuple given in degrees and in ICRS system
    pub fn proj_lonlat(&self, lonlat: &LonLat) -> Option<ImgXY> {
        let lonlat = &self.coo_system.from_icrs(lonlat.clone());

        let img_xy = match &self.proj {
            // Zenithal
            WCSCelestialProj::Azp(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Szp(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Tan(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Stg(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Sin(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Arc(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Zpn(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Zea(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Air(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Ncp(wcs) => wcs.lonlat2img(lonlat),
            // Pseudo-cyl
            WCSCelestialProj::Cyp(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Cea(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Car(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Mer(wcs) => wcs.lonlat2img(lonlat),
            // Cylindrical
            WCSCelestialProj::Sfl(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Par(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Mol(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Ait(wcs) => wcs.lonlat2img(lonlat),
            // Conic
            WCSCelestialProj::Cop(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Cod(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Coe(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::Coo(wcs) => wcs.lonlat2img(lonlat),
            // Hybrid
            WCSCelestialProj::Hpx(wcs) => wcs.lonlat2img(lonlat),

            /* Sip variants */
            // Zenithal
            WCSCelestialProj::AzpSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::SzpSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::TanSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::StgSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::SinSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::ArcSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::ZpnSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::ZeaSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::AirSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::NcpSip(wcs) => wcs.lonlat2img(lonlat),
            // Pseudo-cyl
            WCSCelestialProj::CypSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CeaSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CarSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::MerSip(wcs) => wcs.lonlat2img(lonlat),
            // Cylindrical
            WCSCelestialProj::SflSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::ParSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::MolSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::AitSip(wcs) => wcs.lonlat2img(lonlat),
            // Conic
            WCSCelestialProj::CopSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CodSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CoeSip(wcs) => wcs.lonlat2img(lonlat),
            WCSCelestialProj::CooSip(wcs) => wcs.lonlat2img(lonlat),
            // Hybrid
            WCSCelestialProj::HpxSip(wcs) => wcs.lonlat2img(lonlat),
        };

        img_xy
        //img_xy.map(|xy| ImgXY::new(xy.x() - 1.0, xy.y() - 1.0))
    }

    pub fn proj_xyz(&self, xyz: &(f64, f64, f64)) -> Option<ImgXY> {
        let xyz = &self.coo_system.from_icrs_xyz(XYZ::new(xyz.0, xyz.1, xyz.2));

        let img_xy = match &self.proj {
            // Zenithal
            WCSCelestialProj::Azp(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Szp(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Tan(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Stg(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Sin(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Arc(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Zpn(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Zea(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Air(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Ncp(wcs) => wcs.xyz2img(xyz),
            // Pseudo-cyl
            WCSCelestialProj::Cyp(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Cea(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Car(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Mer(wcs) => wcs.xyz2img(xyz),
            // Cylindrical
            WCSCelestialProj::Sfl(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Par(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Mol(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Ait(wcs) => wcs.xyz2img(xyz),
            // Conic
            WCSCelestialProj::Cop(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Cod(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Coe(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::Coo(wcs) => wcs.xyz2img(xyz),
            // Hybrid
            WCSCelestialProj::Hpx(wcs) => wcs.xyz2img(xyz),

            /* Sip variants */
            // Zenithal
            WCSCelestialProj::AzpSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::SzpSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::TanSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::StgSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::SinSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::ArcSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::ZpnSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::ZeaSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::AirSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::NcpSip(wcs) => wcs.xyz2img(xyz),
            // Pseudo-cyl
            WCSCelestialProj::CypSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CeaSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CarSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::MerSip(wcs) => wcs.xyz2img(xyz),
            // Cylindrical
            WCSCelestialProj::SflSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::ParSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::MolSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::AitSip(wcs) => wcs.xyz2img(xyz),
            // Conic
            WCSCelestialProj::CopSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CodSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CoeSip(wcs) => wcs.xyz2img(xyz),
            WCSCelestialProj::CooSip(wcs) => wcs.xyz2img(xyz),
            // Hybrid
            WCSCelestialProj::HpxSip(wcs) => wcs.xyz2img(xyz),
        };

        img_xy
        //img_xy.map(|xy| ImgXY::new(xy.x() - 1.0, xy.y() - 1.0))
    }

    pub fn unproj_xyz(&self, img_pos: &ImgXY) -> Option<XYZ> {
        //let img_pos = ImgXY::new(img_pos.x() + 1.0, img_pos.y() + 1.0);
        let xyz = match &self.proj {
            // Zenithal
            WCSCelestialProj::Azp(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Szp(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Tan(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Stg(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Sin(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Arc(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Zpn(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Zea(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Air(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Ncp(wcs) => wcs.img2xyz(&img_pos),
            // Pseudo-cyl
            WCSCelestialProj::Cyp(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Cea(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Car(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Mer(wcs) => wcs.img2xyz(&img_pos),
            // Cylindrical
            WCSCelestialProj::Sfl(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Par(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Mol(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Ait(wcs) => wcs.img2xyz(&img_pos),
            // Conic
            WCSCelestialProj::Cop(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Cod(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Coe(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::Coo(wcs) => wcs.img2xyz(&img_pos),
            // Hybrid
            WCSCelestialProj::Hpx(wcs) => wcs.img2xyz(&img_pos),

            /* Sip variants */
            // Zenithal
            WCSCelestialProj::AzpSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::SzpSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::TanSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::StgSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::SinSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::ArcSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::ZpnSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::ZeaSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::AirSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::NcpSip(wcs) => wcs.img2xyz(&img_pos),
            // Pseudo-cyl
            WCSCelestialProj::CypSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::CeaSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::CarSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::MerSip(wcs) => wcs.img2xyz(&img_pos),
            // Cylindrical
            WCSCelestialProj::SflSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::ParSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::MolSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::AitSip(wcs) => wcs.img2xyz(&img_pos),
            // Conic
            WCSCelestialProj::CopSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::CodSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::CoeSip(wcs) => wcs.img2xyz(&img_pos),
            WCSCelestialProj::CooSip(wcs) => wcs.img2xyz(&img_pos),
            // Hybrid
            WCSCelestialProj::HpxSip(wcs) => wcs.img2xyz(&img_pos),
        };
        xyz.map(|v| self.coo_system.to_icrs_xyz(v))
    }

    /// Unproject a (X, Y) point from the image space to get its corresponding location on the sphere
    ///
    /// The result is (lon, lat) tuple expressed in degrees in ICRS
    ///
    /// # Arguments
    ///
    /// * `img_pos`: the image space point expressed as a (X, Y) tuple given en pixels
    pub fn unproj_lonlat(&self, img_pos: &ImgXY) -> Option<LonLat> {
        //let img_pos = ImgXY::new(img_pos.x() + 1.0, img_pos.y() + 1.0);
        let lonlat = match &self.proj {
            // Zenithal
            WCSCelestialProj::Azp(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Szp(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Tan(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Stg(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Sin(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Arc(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Zpn(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Zea(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Air(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Ncp(wcs) => wcs.img2lonlat(&img_pos),
            // Pseudo-cyl
            WCSCelestialProj::Cyp(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Cea(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Car(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Mer(wcs) => wcs.img2lonlat(&img_pos),
            // Cylindrical
            WCSCelestialProj::Sfl(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Par(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Mol(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Ait(wcs) => wcs.img2lonlat(&img_pos),
            // Conic
            WCSCelestialProj::Cop(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Cod(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Coe(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::Coo(wcs) => wcs.img2lonlat(&img_pos),
            // Hybrid
            WCSCelestialProj::Hpx(wcs) => wcs.img2lonlat(&img_pos),

            /* Sip variants */
            // Zenithal
            WCSCelestialProj::AzpSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::SzpSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::TanSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::StgSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::SinSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::ArcSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::ZpnSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::ZeaSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::AirSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::NcpSip(wcs) => wcs.img2lonlat(&img_pos),
            // Pseudo-cyl
            WCSCelestialProj::CypSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::CeaSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::CarSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::MerSip(wcs) => wcs.img2lonlat(&img_pos),
            // Cylindrical
            WCSCelestialProj::SflSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::ParSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::MolSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::AitSip(wcs) => wcs.img2lonlat(&img_pos),
            // Conic
            WCSCelestialProj::CopSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::CodSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::CoeSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSCelestialProj::CooSip(wcs) => wcs.img2lonlat(&img_pos),
            // Hybrid
            WCSCelestialProj::HpxSip(wcs) => wcs.img2lonlat(&img_pos),
        };

        lonlat.map(|ll| self.coo_system.to_icrs(ll))
    }

    /// Getter of the coordinate system
    pub fn coo_system(&self) -> &CooSystem {
        &self.coo_system
    }
}

#[cfg(test)]
mod tests {
    use super::WCS;
    use crate::WCSParams;
    use crate::Error;

    use crate::mapproj::Projection;
    use fitsrs::fits::Fits;
    use fitsrs::card::CardValue;
    use fitsrs::hdu::header::{extension::image::Image, Header};
    use fitsrs::ImageData;
    use fitsrs::card::Value;

    use glob::glob;
    use mapproj::{CanonicalProjection, ImgXY, LonLat};
    use std::f64::consts::PI;
    use std::fs::File;
    use std::io::BufReader;

    use std::convert::TryFrom;
    use std::str::FromStr;
    
    fn parse_optional_card_with_type<T: CardValue + FromStr>(header: &Header<Image>, key: &'static str) -> Result<Option<T>, Error> {
        match  header.get_parsed::<T>(key).transpose() {
            Ok(v) => Ok(v),
            _ => {
                let str = header.get_parsed::<String>(key).transpose()
                    .unwrap_or(None);
    
                Ok(if let Some(ss) = str {
                    ss.trim().parse::<T>()
                        .map(|v| Some(v))
                        .unwrap_or(None)
                } else {
                    // card not found but it is ok as it is not mandatory
                    None
                })
            }
        }
    }
    
    fn parse_mandatory_card_with_type<T: CardValue>(header: &Header<Image>, key: &'static str) -> Result<T, &'static str> {
        match header.get_parsed::<T>(key) {
            // No parsing error and found
            Some(Ok(v)) => {
                Ok(v)
            },
            // No error but not found, we return an error
            None => Err("Mandatory keyword not found"),
            // Return the parsing error
            Some(Err(_)) => {
                Err("Error parsing mandatory keyword")
            }
        }
    }
    
    impl<'a> TryFrom<&'a Header<Image>> for WCS {
        type Error = Error;
    
        fn try_from(h: &'a Header<Image>) -> Result<Self, Self::Error> {
            let params = WCSParams {
                naxis: parse_mandatory_card_with_type::<i64>(h, "NAXIS").unwrap(),
                ctype1: parse_mandatory_card_with_type::<String>(h, "CTYPE1").unwrap(),

                naxis1: parse_optional_card_with_type::<i64>(h, "NAXIS1")?,
                naxis2: parse_optional_card_with_type::<i64>(h, "NAXIS2")?,

                ctype2: parse_optional_card_with_type::<String>(h, "CTYPE2")?,
                ctype3: parse_optional_card_with_type::<String>(h, "CTYPE3")?,

                a_order: parse_optional_card_with_type::<i64>(h, "A_ORDER")?,
                b_order: parse_optional_card_with_type::<i64>(h, "B_ORDER")?,
                ap_order: parse_optional_card_with_type::<i64>(h, "AP_ORDER")?,
                bp_order: parse_optional_card_with_type::<i64>(h, "BP_ORDER")?,
                crpix1: parse_optional_card_with_type::<f64>(h, "CRPIX1")?,
                crpix2: parse_optional_card_with_type::<f64>(h, "CRPIX2")?,
                crpix3: parse_optional_card_with_type::<f64>(h, "CRPIX3")?,
                crval1: parse_optional_card_with_type::<f64>(h, "CRVAL1")?,
                crval2: parse_optional_card_with_type::<f64>(h, "CRVAL2")?,
                crval3: parse_optional_card_with_type::<f64>(h, "CRVAL3")?,
                crota1: parse_optional_card_with_type::<f64>(h, "CROTA1")?,
                crota2: parse_optional_card_with_type::<f64>(h, "CROTA2")?,
                crota3: parse_optional_card_with_type::<f64>(h, "CROTA3")?,
                cdelt1: parse_optional_card_with_type::<f64>(h, "CDELT1")?,
                cdelt2: parse_optional_card_with_type::<f64>(h, "CDELT2")?,
                cdelt3: parse_optional_card_with_type::<f64>(h, "CDELT3")?,
                naxis3: parse_optional_card_with_type::<i64>(h, "NAXIS3")?,
                naxis4: parse_optional_card_with_type::<i64>(h, "NAXIS4")?,
                lonpole: parse_optional_card_with_type::<f64>(h, "LONPOLE")?,
                latpole: parse_optional_card_with_type::<f64>(h, "LATPOLE")?,
                equinox: parse_optional_card_with_type::<f64>(h, "EQUINOX")?,
                epoch: parse_optional_card_with_type::<f64>(h, "EPOCH")?,
                radesys: parse_optional_card_with_type::<String>(h, "RADESYS")?,
                pv1_0: parse_optional_card_with_type::<f64>(h, "PV1_0")?,
                pv1_1: parse_optional_card_with_type::<f64>(h, "PV1_1")?,
                pv1_2: parse_optional_card_with_type::<f64>(h, "PV1_2")?,
                pv2_0: parse_optional_card_with_type::<f64>(h, "PV2_0")?,
                pv2_1: parse_optional_card_with_type::<f64>(h, "PV2_1")?,
                pv2_2: parse_optional_card_with_type::<f64>(h, "PV2_2")?,
                pv2_3: parse_optional_card_with_type::<f64>(h, "PV2_3")?,
                pv2_4: parse_optional_card_with_type::<f64>(h, "PV2_4")?,
                pv2_5: parse_optional_card_with_type::<f64>(h, "PV2_5")?,
                pv2_6: parse_optional_card_with_type::<f64>(h, "PV2_6")?,
                pv2_7: parse_optional_card_with_type::<f64>(h, "PV2_7")?,
                pv2_8: parse_optional_card_with_type::<f64>(h, "PV2_8")?,
                pv2_9: parse_optional_card_with_type::<f64>(h, "PV2_9")?,
                pv2_10: parse_optional_card_with_type::<f64>(h, "PV2_10")?,
                pv2_11: parse_optional_card_with_type::<f64>(h, "PV2_11")?,
                pv2_12: parse_optional_card_with_type::<f64>(h, "PV2_12")?,
                pv2_13: parse_optional_card_with_type::<f64>(h, "PV2_13")?,
                pv2_14: parse_optional_card_with_type::<f64>(h, "PV2_14")?,
                pv2_15: parse_optional_card_with_type::<f64>(h, "PV2_15")?,
                pv2_16: parse_optional_card_with_type::<f64>(h, "PV2_16")?,
                pv2_17: parse_optional_card_with_type::<f64>(h, "PV2_17")?,
                pv2_18: parse_optional_card_with_type::<f64>(h, "PV2_18")?,
                pv2_19: parse_optional_card_with_type::<f64>(h, "PV2_19")?,
                pv2_20: parse_optional_card_with_type::<f64>(h, "PV2_20")?,
                cd1_1: parse_optional_card_with_type::<f64>(h, "CD1_1")?,
                cd1_2: parse_optional_card_with_type::<f64>(h, "CD1_2")?,
                cd1_3: parse_optional_card_with_type::<f64>(h, "CD1_3")?,
                cd2_1: parse_optional_card_with_type::<f64>(h, "CD2_1")?,
                cd2_2: parse_optional_card_with_type::<f64>(h, "CD2_2")?,
                cd2_3: parse_optional_card_with_type::<f64>(h, "CD2_3")?,
                cd3_1: parse_optional_card_with_type::<f64>(h, "CD3_1")?,
                cd3_2: parse_optional_card_with_type::<f64>(h, "CD3_2")?,
                cd3_3: parse_optional_card_with_type::<f64>(h, "CD3_3")?,
                pc1_1: parse_optional_card_with_type::<f64>(h, "PC1_1")?,
                pc1_2: parse_optional_card_with_type::<f64>(h, "PC1_2")?,
                pc1_3: parse_optional_card_with_type::<f64>(h, "PC1_3")?,
                pc2_1: parse_optional_card_with_type::<f64>(h, "PC2_1")?,
                pc2_2: parse_optional_card_with_type::<f64>(h, "PC2_2")?,
                pc2_3: parse_optional_card_with_type::<f64>(h, "PC2_3")?,
                pc3_1: parse_optional_card_with_type::<f64>(h, "PC3_1")?,
                pc3_2: parse_optional_card_with_type::<f64>(h, "PC3_2")?,
                pc3_3: parse_optional_card_with_type::<f64>(h, "PC3_3")?,
                a_0_0: parse_optional_card_with_type::<f64>(h, "A_0_0")?,
                a_1_0: parse_optional_card_with_type::<f64>(h, "A_1_0")?,
                a_2_0: parse_optional_card_with_type::<f64>(h, "A_2_0")?,
                a_3_0: parse_optional_card_with_type::<f64>(h, "A_3_0")?,
                a_4_0: parse_optional_card_with_type::<f64>(h, "A_4_0")?,
                a_5_0: parse_optional_card_with_type::<f64>(h, "A_5_0")?,
                a_6_0: parse_optional_card_with_type::<f64>(h, "A_6_0")?,
                a_0_1: parse_optional_card_with_type::<f64>(h, "A_0_1")?,
                a_1_1: parse_optional_card_with_type::<f64>(h, "A_1_1")?,
                a_2_1: parse_optional_card_with_type::<f64>(h, "A_2_1")?,
                a_3_1: parse_optional_card_with_type::<f64>(h, "A_3_1")?,
                a_4_1: parse_optional_card_with_type::<f64>(h, "A_4_1")?,
                a_5_1: parse_optional_card_with_type::<f64>(h, "A_5_1")?,
                a_0_2: parse_optional_card_with_type::<f64>(h, "A_0_2")?,
                a_1_2: parse_optional_card_with_type::<f64>(h, "A_1_2")?,
                a_2_2: parse_optional_card_with_type::<f64>(h, "A_2_2")?,
                a_3_2: parse_optional_card_with_type::<f64>(h, "A_3_2")?,
                a_4_2: parse_optional_card_with_type::<f64>(h, "A_4_2")?,
                a_0_3: parse_optional_card_with_type::<f64>(h, "A_0_3")?,
                a_1_3: parse_optional_card_with_type::<f64>(h, "A_1_3")?,
                a_2_3: parse_optional_card_with_type::<f64>(h, "A_2_3")?,
                a_3_3: parse_optional_card_with_type::<f64>(h, "A_3_3")?,
                a_0_4: parse_optional_card_with_type::<f64>(h, "A_0_4")?,
                a_1_4: parse_optional_card_with_type::<f64>(h, "A_1_4")?,
                a_2_4: parse_optional_card_with_type::<f64>(h, "A_2_4")?,
                a_0_5: parse_optional_card_with_type::<f64>(h, "A_0_5")?,
                a_1_5: parse_optional_card_with_type::<f64>(h, "A_1_5")?,
                a_0_6: parse_optional_card_with_type::<f64>(h, "A_0_6")?,
                ap_0_0: parse_optional_card_with_type::<f64>(h, "AP_0_0")?,
                ap_1_0: parse_optional_card_with_type::<f64>(h, "AP_1_0")?,
                ap_2_0: parse_optional_card_with_type::<f64>(h, "AP_2_0")?,
                ap_3_0: parse_optional_card_with_type::<f64>(h, "AP_3_0")?,
                ap_4_0: parse_optional_card_with_type::<f64>(h, "AP_4_0")?,
                ap_5_0: parse_optional_card_with_type::<f64>(h, "AP_5_0")?,
                ap_6_0: parse_optional_card_with_type::<f64>(h, "AP_6_0")?,
                ap_0_1: parse_optional_card_with_type::<f64>(h, "AP_0_1")?,
                ap_1_1: parse_optional_card_with_type::<f64>(h, "AP_1_1")?,
                ap_2_1: parse_optional_card_with_type::<f64>(h, "AP_2_1")?,
                ap_3_1: parse_optional_card_with_type::<f64>(h, "AP_3_1")?,
                ap_4_1: parse_optional_card_with_type::<f64>(h, "AP_4_1")?,
                ap_5_1: parse_optional_card_with_type::<f64>(h, "AP_5_1")?,
                ap_0_2: parse_optional_card_with_type::<f64>(h, "AP_0_2")?,
                ap_1_2: parse_optional_card_with_type::<f64>(h, "AP_1_2")?,
                ap_2_2: parse_optional_card_with_type::<f64>(h, "AP_2_2")?,
                ap_3_2: parse_optional_card_with_type::<f64>(h, "AP_3_2")?,
                ap_4_2: parse_optional_card_with_type::<f64>(h, "AP_4_2")?,
                ap_0_3: parse_optional_card_with_type::<f64>(h, "AP_0_3")?,
                ap_1_3: parse_optional_card_with_type::<f64>(h, "AP_1_3")?,
                ap_2_3: parse_optional_card_with_type::<f64>(h, "AP_2_3")?,
                ap_3_3: parse_optional_card_with_type::<f64>(h, "AP_3_3")?,
                ap_0_4: parse_optional_card_with_type::<f64>(h, "AP_0_4")?,
                ap_1_4: parse_optional_card_with_type::<f64>(h, "AP_1_4")?,
                ap_2_4: parse_optional_card_with_type::<f64>(h, "AP_2_4")?,
                ap_0_5: parse_optional_card_with_type::<f64>(h, "AP_0_5")?,
                ap_1_5: parse_optional_card_with_type::<f64>(h, "AP_1_5")?,
                ap_0_6: parse_optional_card_with_type::<f64>(h, "AP_0_6")?,
                b_0_0: parse_optional_card_with_type::<f64>(h, "B_0_0")?,
                b_1_0: parse_optional_card_with_type::<f64>(h, "B_1_0")?,
                b_2_0: parse_optional_card_with_type::<f64>(h, "B_2_0")?,
                b_3_0: parse_optional_card_with_type::<f64>(h, "B_3_0")?,
                b_4_0: parse_optional_card_with_type::<f64>(h, "B_4_0")?,
                b_5_0: parse_optional_card_with_type::<f64>(h, "B_5_0")?,
                b_6_0: parse_optional_card_with_type::<f64>(h, "B_6_0")?,
                b_0_1: parse_optional_card_with_type::<f64>(h, "B_0_1")?,
                b_1_1: parse_optional_card_with_type::<f64>(h, "B_1_1")?,
                b_2_1: parse_optional_card_with_type::<f64>(h, "B_2_1")?,
                b_3_1: parse_optional_card_with_type::<f64>(h, "B_3_1")?,
                b_4_1: parse_optional_card_with_type::<f64>(h, "B_4_1")?,
                b_5_1: parse_optional_card_with_type::<f64>(h, "B_5_1")?,
                b_0_2: parse_optional_card_with_type::<f64>(h, "B_0_2")?,
                b_1_2: parse_optional_card_with_type::<f64>(h, "B_1_2")?,
                b_2_2: parse_optional_card_with_type::<f64>(h, "B_2_2")?,
                b_3_2: parse_optional_card_with_type::<f64>(h, "B_3_2")?,
                b_4_2: parse_optional_card_with_type::<f64>(h, "B_4_2")?,
                b_0_3: parse_optional_card_with_type::<f64>(h, "B_0_3")?,
                b_1_3: parse_optional_card_with_type::<f64>(h, "B_1_3")?,
                b_2_3: parse_optional_card_with_type::<f64>(h, "B_2_3")?,
                b_3_3: parse_optional_card_with_type::<f64>(h, "B_3_3")?,
                b_0_4: parse_optional_card_with_type::<f64>(h, "B_0_4")?,
                b_1_4: parse_optional_card_with_type::<f64>(h, "B_1_4")?,
                b_2_4: parse_optional_card_with_type::<f64>(h, "B_2_4")?,
                b_0_5: parse_optional_card_with_type::<f64>(h, "B_0_5")?,
                b_1_5: parse_optional_card_with_type::<f64>(h, "B_1_5")?,
                b_0_6: parse_optional_card_with_type::<f64>(h, "B_0_6")?,
                bp_0_0: parse_optional_card_with_type::<f64>(h, "BP_0_0")?,
                bp_1_0: parse_optional_card_with_type::<f64>(h, "BP_1_0")?,
                bp_2_0: parse_optional_card_with_type::<f64>(h, "BP_2_0")?,
                bp_3_0: parse_optional_card_with_type::<f64>(h, "BP_3_0")?,
                bp_4_0: parse_optional_card_with_type::<f64>(h, "BP_4_0")?,
                bp_5_0: parse_optional_card_with_type::<f64>(h, "BP_5_0")?,
                bp_6_0: parse_optional_card_with_type::<f64>(h, "BP_6_0")?,
                bp_0_1: parse_optional_card_with_type::<f64>(h, "BP_0_1")?,
                bp_1_1: parse_optional_card_with_type::<f64>(h, "BP_1_1")?,
                bp_2_1: parse_optional_card_with_type::<f64>(h, "BP_2_1")?,
                bp_3_1: parse_optional_card_with_type::<f64>(h, "BP_3_1")?,
                bp_4_1: parse_optional_card_with_type::<f64>(h, "BP_4_1")?,
                bp_5_1: parse_optional_card_with_type::<f64>(h, "BP_5_1")?,
                bp_0_2: parse_optional_card_with_type::<f64>(h, "BP_0_2")?,
                bp_1_2: parse_optional_card_with_type::<f64>(h, "BP_1_2")?,
                bp_2_2: parse_optional_card_with_type::<f64>(h, "BP_2_2")?,
                bp_3_2: parse_optional_card_with_type::<f64>(h, "BP_3_2")?,
                bp_4_2: parse_optional_card_with_type::<f64>(h, "BP_4_2")?,
                bp_0_3: parse_optional_card_with_type::<f64>(h, "BP_0_3")?,
                bp_1_3: parse_optional_card_with_type::<f64>(h, "BP_1_3")?,
                bp_2_3: parse_optional_card_with_type::<f64>(h, "BP_2_3")?,
                bp_3_3: parse_optional_card_with_type::<f64>(h, "BP_3_3")?,
                bp_0_4: parse_optional_card_with_type::<f64>(h, "BP_0_4")?,
                bp_1_4: parse_optional_card_with_type::<f64>(h, "BP_1_4")?,
                bp_2_4: parse_optional_card_with_type::<f64>(h, "BP_2_4")?,
                bp_0_5: parse_optional_card_with_type::<f64>(h, "BP_0_5")?,
                bp_1_5: parse_optional_card_with_type::<f64>(h, "BP_1_5")?,
                bp_0_6: parse_optional_card_with_type::<f64>(h, "BP_0_6")?,
            };
    
            WCS::new(&params)
        }
    }
    
    fn wcs_from_fits_header(header: &Header<Image>) -> Result<WCS, Error> {
        header.try_into()
    }

    #[test]
    fn test_visualize() {
        let f = File::open("examples/panstarrs-rotated-around-orion.fits").unwrap();

        let reader = BufReader::new(f);
        let mut fits = Fits::from_reader(reader);
        let hdu = fits.next().unwrap().unwrap();

        match hdu {
            HDU::XImage(hdu) | HDU::Primary(hdu) => {
                let header = hdu.get_header();

                // Parse data
                let data = match fits.get_data(&hdu) {
                    ImageData::F32(it) => it.collect::<Vec<_>>(),
                    _ => unreachable!(),
                };

                let wcs = wcs_from_fits_header(&header).unwrap();
                reproject_fits_image(mapproj::zenithal::azp::Azp::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::zenithal::szp::Szp::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::zenithal::tan::Tan::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::zenithal::stg::Stg::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::zenithal::sin::Sin::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::zenithal::arc::Arc::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::zenithal::zea::Zea::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::zenithal::air::Air::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::zenithal::ncp::Ncp::new(), &wcs, &header, &data);

                reproject_fits_image(mapproj::pseudocyl::mol::Mol::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::pseudocyl::ait::Ait::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::pseudocyl::par::Par::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::pseudocyl::sfl::Sfl::new(), &wcs, &header, &data);

                reproject_fits_image(mapproj::cylindrical::cyp::Cyp::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::cylindrical::cea::Cea::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::cylindrical::car::Car::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::cylindrical::mer::Mer::new(), &wcs, &header, &data);

                reproject_fits_image(mapproj::conic::cod::Cod::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::conic::cop::Cop::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::conic::coo::Coo::new(), &wcs, &header, &data);
                reproject_fits_image(mapproj::conic::coe::Coe::new(), &wcs, &header, &data);

                reproject_fits_image(mapproj::hybrid::hpx::Hpx::new(), &wcs, &header, &data);
            }
            _ => unreachable!()
        }
    }

    fn reproject_fits_image<'a, T: CanonicalProjection>(
        proj: T,
        wcs: &WCS,
        header: &Header<Image>,
        data: &[f32],
    ) {
        let scale = header
            .get_parsed::<f64>("BSCALE")
            .unwrap_or(Ok(1.0))
            .unwrap() as f32;
        let offset = header
            .get_parsed::<f64>("BZERO")
            .unwrap_or(Ok(0.0))
            .unwrap() as f32;

        let xtension = header.get_xtension();
        let width = *xtension.get_naxisn(1).unwrap();
        let height = *xtension.get_naxisn(2).unwrap();

        //let proj = mapproj::zenithal::sin::Sin::new();
        let bounds = proj.bounds();
        let x_bounds = bounds.x_bounds().as_ref().unwrap_or(&((-PI)..=PI));
        let y_bounds = bounds.y_bounds().as_ref().unwrap_or(&((-PI)..=PI));

        let x_off = x_bounds.start();
        let x_len = x_bounds.end() - x_bounds.start();

        let y_off = y_bounds.start();
        let y_len = y_bounds.end() - y_bounds.start();

        const WIDTH_IMAGE: usize = 1024;
        const HEIGHT_IMAGE: usize = 1024;
        // Create a new ImgBuf with width: imgx and height: imgy
        let mut imgbuf = image::ImageBuffer::new(WIDTH_IMAGE as u32, HEIGHT_IMAGE as u32);

        for y in 0..height {
            for x in 0..width {
                let grayscale_val = (data[(y * width + x) as usize] * scale + offset) as u8;

                let img_xy = ImgXY::new(x as f64, y as f64);
                if let Some(lonlat) = wcs.unproj(&img_xy) {
                    if let Some(proj_xy) = proj.proj_lonlat(&lonlat) {
                        let proj_x = ((proj_xy.x() as f64) - x_off) / x_len; // between 0 and 1
                        let proj_y = ((proj_xy.y() as f64) - y_off) / y_len; // between 0 and 1

                        if (0.0..1.0).contains(&proj_x) && (0.0..1.0).contains(&proj_y) {
                            let ix = (proj_x * (WIDTH_IMAGE as f64)) as usize;
                            let iy = (proj_y * (HEIGHT_IMAGE as f64)) as usize;

                            let pixel = imgbuf
                                .get_pixel_mut(ix as u32, (HEIGHT_IMAGE as u32) - iy as u32 - 1);
                            *pixel = image::Rgb([grayscale_val, grayscale_val, grayscale_val]);
                        }
                    }
                }
            }
        }

        let filename = &format!(
            "tests/reproj/pans-{}.jpeg",
            <T as CanonicalProjection>::WCS_NAME
        );
        imgbuf.save(filename).unwrap();
    }

    macro_rules! assert_delta {
        ($x:expr, $y:expr, $d:expr) => {
            if ($x - $y).abs() > $d {
                panic!();
            }
        };
    }

    #[test]
    fn astropy_comparison() {
        use std::fs;

        let fits_file_paths = fs::read_dir("./examples")
            .unwrap()
            .filter_map(|res| res.ok())
            // Map the directory entries to paths
            .map(|dir_entry| dir_entry.path())
            // Filter out all paths with extensions other than `csv`
            .filter_map(|path| {
                if path.extension().map_or(false, |ext| ext == "fits") {
                    Some(path)
                } else {
                    None
                }
            });

        for path in fits_file_paths {
            println!("Test {:?}", path.display());

            let f = File::open(path.clone()).unwrap();

            let reader = BufReader::new(f);
            let mut fits = Fits::from_reader(reader);
            let hdu = fits.next().unwrap().unwrap();
            match hdu {
                HDU::XImage(hdu) | HDU::Primary(hdu) => {
                    let header = hdu.get_header();
                    let wcs = wcs_from_fits_header(header).unwrap();
        
                    // add the astropy opened results
                    let path_astropy = path.with_extension("fits.csv");
                    let f = File::open(path_astropy).unwrap();
        
                    let mut rdr = csv::Reader::from_reader(BufReader::new(f));
                    for result in rdr.records() {
                        let record = result.unwrap();
        
                        let ra: f64 = record[0].parse().unwrap();
                        let dec: f64 = record[1].parse().unwrap();
                        let x: f64 = record[2].parse().unwrap();
                        let y: f64 = record[3].parse().unwrap();
        
                        if ra.is_finite() && dec.is_finite() {
                            if let Some(img_xy) = wcs.proj(&LonLat::new(ra, dec)) {
                                //dbg!(img_xy.x() - x);
                                //dbg!(img_xy.y() - y);
        
                                assert_delta!(img_xy.x(), x, 1e-4);
                                assert_delta!(img_xy.y(), y, 1e-4);
                            }
                        }
                    }
                }
                _ => unreachable!()
            };
        }
    }

    #[test]
    fn crval_to_crpix() {
        for entry in glob("examples/*.fits").unwrap() {
            if let Ok(path) = dbg!(entry) {
                let f = File::open(path).unwrap();
                let reader = BufReader::new(f);
                let mut fits = Fits::from_reader(reader);
                let hdu = fits.next().unwrap().unwrap();

                match hdu {
                    HDU::XImage(hdu) | HDU::Primary(hdu) => {
                        let header = hdu.get_header();
                        let crval1 = header
                            .get_parsed::<f64>("CRVAL1")
                            .unwrap_or(Ok(0.0))
                            .unwrap();
                        let crval2 = header
                            .get_parsed::<f64>("CRVAL2")
                            .unwrap_or(Ok(0.0))
                            .unwrap();
                        let crpix1 = if let Some(Value::Integer { value, .. }) = header.get("CRPIX1") {
                            *value as f64
                        } else if let Some(Value::Float { value, .. }) = header.get("CRPIX1") {
                            *value
                        } else {
                            0.0
                        };
                        
                        let crpix2 = if let Some(Value::Integer { value, .. }) = header.get("CRPIX2") {
                            *value as f64
                        } else if let Some(Value::Float { value, .. }) = header.get("CRPIX2") {
                            *value
                        } else {
                            0.0
                        };
        
                        let wcs = wcs_from_fits_header(&header).unwrap();
        
                        // crval to crpix
                        let proj_px = wcs
                            .proj(&LonLat::new(
                                dbg!(crval1).to_radians(),
                                dbg!(crval2).to_radians(),
                            ))
                            .unwrap();
                        assert_delta!(proj_px.x(), crpix1, 1e-6);
                        assert_delta!(proj_px.y(), crpix2, 1e-6);
        
                        // crpix to crval
                        let lonlat = wcs
                            .unproj_lonlat(&ImgXY::new(dbg!(crpix1), dbg!(crpix2)))
                            .unwrap();
                        assert_delta!(lonlat.lon(), crval1.to_radians(), 1e-6);
                        assert_delta!(lonlat.lat(), crval2.to_radians(), 1e-6);
                    },
                    _ => unreachable!()
                };
            }
        }
    }

    use fitsrs::hdu::HDU;
    #[test]
    fn open_fits() {
        let f = File::open("examples/neowise.fits").unwrap();

        let reader = BufReader::new(f);
        let mut fits = Fits::from_reader(reader);
        let hdu = fits.next().unwrap().unwrap();
        match hdu {
            HDU::XImage(hdu) | HDU::Primary(hdu) => {
                let header = hdu.get_header();
                let wcs = wcs_from_fits_header(header).unwrap();
                dbg!(wcs.unproj(&ImgXY::new(0.0, 1200.0)));
            },
            _ => unreachable!()
        }
    }
}
