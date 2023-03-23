extern crate fitsrs;
extern crate mapproj;
#[macro_use]
extern crate quick_error;

mod error;

use error::Error;
mod projection;
mod sip;
mod utils;

use crate::projection::WCSCanonicalProjection;

// Imports
use fitsrs::hdu::header::{
    Header,
    extension::image::Image
};
use mapproj::{
    conic::{cod::Cod, coe::Coe, coo::Coo, cop::Cop},
    cylindrical::{car::Car, cea::Cea, cyp::Cyp, mer::Mer},
    img2celestial::Img2Celestial,
    img2proj::{WcsImgXY2ProjXY, WcsWithSipImgXY2ProjXY},
    pseudocyl::{ait::Ait, mol::Mol, par::Par, sfl::Sfl},
    zenithal::{
        air::Air, arc::Arc, azp::Azp, sin::Sin, stg::Stg, szp::Szp, tan::Tan, zea::Zea, zpn::Zpn, ncp::Ncp
    },
};

use paste::paste;
/// macro
macro_rules! create_specific_proj {
    ( $proj_name:ident, $header:expr, $ctype1:expr, $crpix1:expr, $crpix2:expr, $img2proj:expr ) => {{
        let proj = $proj_name::parse_proj(&$header)?;

        let is_sip_found = &$ctype1[($ctype1.len() - 3)..] == "SIP";
        if is_sip_found {
            let sip = sip::parse_sip($header, $crpix1, $crpix2)?;
            let img2proj = WcsWithSipImgXY2ProjXY::new($img2proj, sip);

            paste! {
                Ok(WCSProj::[ <$proj_name Sip> ](Img2Celestial::new(img2proj, proj)))
            }
        } else {
            Ok(WCSProj::$proj_name(Img2Celestial::new($img2proj, proj)))
        }
    }};
}

/// Structure alias coming from mapproj defining
/// image space pixel coordinates
pub type ImgXY = mapproj::ImgXY;
/// Structure alias coming from mapproj defining
/// longitude and latitude expressed in degrees
pub type LonLat = mapproj::LonLat;


pub struct WCS {
    /* Metadata keywords */
    /// Width of the image in pixels
    naxis1: u64,
    /// Height of the image in pixels
    naxis2: u64,
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
    /// Create a WCS from a specific fits header parsed with fitsrs
    /// # Param
    /// * `header`: Header unit coming from fitsrs.
    ///   This contains all the cards of one HDU.
    pub fn new(header: &Header<Image>) -> Result<Self, Error> {
        let xtension = header.get_xtension();

        let naxis1 = xtension.get_naxisn(1).ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS1"))?;
        let naxis2 = xtension.get_naxisn(2).ok_or(Error::MandatoryWCSKeywordsMissing("NAXIS2"))?;

        let proj = WCSProj::new(header)?;

        // Compute the field of view along the naxis1 and naxis2 axis
        let center = proj.unproj_lonlat(&ImgXY::new((*naxis1 as f64) / 2.0, (*naxis2 as f64) / 2.0))
            .ok_or(Error::UnprojNotDefined((*naxis1 as f64) / 2.0, (*naxis2 as f64) / 2.0))?;

        let half_fov1 = if let Some(top) = proj.unproj_lonlat(&ImgXY::new((*naxis1 as f64) / 2.0, *naxis2 as f64)) {
            utils::angular_dist(
                top.into(),
                center.clone().into()
            )
        } else {
            180.0_f64.to_radians()
        };

        let half_fov2 = if let Some(left) = proj.unproj_lonlat(&ImgXY::new(0.0, (*naxis2 as f64) / 2.0)) {
            utils::angular_dist(
                left.into(),
                center.into()
            )
        } else {
            180.0_f64.to_radians()
        };

        Ok(WCS {
            naxis1: *naxis1 as u64,
            naxis2: *naxis2 as u64,
            fov1: half_fov1 * 2.0,
            fov2: half_fov2 * 2.0,
            proj: proj,
        })
    }

    /// Returns the dimensions of the image given by the NAXIS1 x NAXIS2 keyword
    pub fn img_dimensions(&self) -> (u64, u64) {
        (self.naxis1, self.naxis2)
    }

    pub fn field_of_view(&self) -> (f64, f64) {
        (self.fov1, self.fov2)
    }

    /// Project a (lon, lat) 3D sphere position to get its corresponding location on the image
    /// The result is given a (X, Y) tuple expressed in pixel coordinates.
    /// 
    /// # Param
    /// * `lonlat`: the 3D sphere vertex expressed as a (lon, lat) tuple given in degrees
    pub fn proj(&self, lonlat: &LonLat) -> Option<ImgXY> {
        self.proj.proj_lonlat(lonlat)
    }

    /// Unproject a (X, Y) point from the image space to get its corresponding location on the sphere
    /// The result is given a (lon, lat) tuple expressed in degrees.
    /// 
    /// # Param
    /// * `img_pos`: the image space point expressed as a (X, Y) tuple given en pixels
    pub fn unproj(&self, img_pos: &ImgXY) -> Option<LonLat> {
        self.proj.unproj_lonlat(img_pos)
    }
}

use std::ops::Deref;
impl Deref for WCS {
    type Target = WCSProj;

    fn deref(&self) -> &Self::Target {
        &self.proj
    }
}

/// Main enum structure descripting a WCS object
/// Once created, the user can proceed two operation on it
/// * The projection of a (lon, lat) tuple onto the image space.
///   Results are given in pixels
/// * The unprojection of a (x, y) tuple given in pixel coordinates onto the sphere.
///   Results are given as a (lon, lat) tuple expressed in degrees
pub enum WCSProj {
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
}

impl WCSProj {
    /// Create a WCS from a specific fits header parsed with fitsrs
    /// # Param
    /// * `header`: Header unit coming from fitsrs.
    ///   This contains all the cards of one HDU.
    pub fn new(header: &Header<Image>) -> Result<Self, Error> {
        // 1. Identify the image <-> intermediate projection
        // a. Linear transformation matrix cases:
        // - CRPIXi + CDij
        // - CRPIXi + CDELTi + CROTA2
        // - CRPIXi + CDELTi + PCij
        let crpix1 = header.get_parsed::<f64>(b"CRPIX1  ").unwrap_or(Ok(0.0))?;
        let crpix2 = header.get_parsed::<f64>(b"CRPIX2  ").unwrap_or(Ok(0.0))?;

        let cdelt1 = header.get_parsed::<f64>(b"CDELT1  ");

        let img2proj = if let Some(cdelt1) = cdelt1 {
            let cdelt1 = cdelt1?;
            let cdelt2 = header.get_parsed::<f64>(b"CDELT2  ").unwrap_or(Ok(1.0))?;

            let pc11 = header.get_parsed::<f64>(b"PC1_1   ");

            if let Some(pc11) = pc11 {
                let pc11 = pc11?;
                // CDELTi + PCij case
                let pc12 = header.get_parsed::<f64>(b"PC1_2   ").unwrap_or(Ok(0.0))?;
                let pc21 = header.get_parsed::<f64>(b"PC2_1   ").unwrap_or(Ok(0.0))?;
                let pc22 = header.get_parsed::<f64>(b"PC2_2   ").unwrap_or(Ok(1.0))?;

                WcsImgXY2ProjXY::from_pc(crpix1, crpix2, pc11, pc12, pc21, pc22, cdelt1, cdelt2)
            } else {
                // CDELTi + CROTA2 case
                let crota2 = header.get_parsed::<f64>(b"CROTA2  ").unwrap_or(Ok(0.0))?;

                WcsImgXY2ProjXY::from_cr(crpix1, crpix2, crota2, cdelt1, cdelt2)
            }
        } else {
            // CDij case
            let cd11 = header.get_parsed::<f64>(b"CD1_1   ").unwrap_or(Ok(1.0))?;
            let cd12 = header.get_parsed::<f64>(b"CD1_2   ").unwrap_or(Ok(0.0))?;
            let cd21 = header.get_parsed::<f64>(b"CD2_1   ").unwrap_or(Ok(0.0))?;
            let cd22 = header.get_parsed::<f64>(b"CD2_2   ").unwrap_or(Ok(1.0))?;

            WcsImgXY2ProjXY::from_cd(crpix1, crpix2, cd11, cd12, cd21, cd22)
        };

        // 2. Identify the projection type
        let ctype1 = utils::retrieve_mandatory_parsed_keyword::<String>(header, "CTYPE1  ")?;
        let _ = utils::retrieve_mandatory_parsed_keyword::<String>(header, "CTYPE2  ")?;

        let proj_name = &ctype1[5..=7];

        match proj_name.as_bytes() {
            // Zenithal
            b"AZP" => {
                create_specific_proj!(Azp, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"SZP" => {
                create_specific_proj!(Szp, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"TAN" => {
                create_specific_proj!(Tan, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"STG" => {
                create_specific_proj!(Stg, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"SIN" => {
                create_specific_proj!(Sin, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"ARC" => {
                create_specific_proj!(Arc, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"ZPN" => {
                create_specific_proj!(Zpn, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"ZEA" => {
                create_specific_proj!(Zea, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"AIR" => {
                create_specific_proj!(Air, header, ctype1, crpix1, crpix2, img2proj)
            }
            // Cylindrical
            b"CYP" => {
                create_specific_proj!(Cyp, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"CEA" => {
                create_specific_proj!(Cea, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"CAR" => {
                create_specific_proj!(Car, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"MER" => {
                create_specific_proj!(Mer, header, ctype1, crpix1, crpix2, img2proj)
            }
            // Pseudo-cylindrical
            b"SFL" => {
                create_specific_proj!(Sfl, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"PAR" => {
                create_specific_proj!(Par, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"MOL" => {
                create_specific_proj!(Mol, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"AIT" => {
                create_specific_proj!(Ait, header, ctype1, crpix1, crpix2, img2proj)
            }
            // Conic
            b"COP" => {
                create_specific_proj!(Cop, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"COD" => {
                create_specific_proj!(Cod, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"COE" => {
                create_specific_proj!(Coe, header, ctype1, crpix1, crpix2, img2proj)
            }
            b"COO" => {
                create_specific_proj!(Coo, header, ctype1, crpix1, crpix2, img2proj)
            }
            _ => Err(Error::NotImplementedProjection(proj_name.to_string())),
        }
    }

    /// Project a (lon, lat) 3D sphere position to get its corresponding location on the image
    /// The result is given a (X, Y) tuple expressed in pixel coordinates.
    /// 
    /// # Param
    /// * `lonlat`: the 3D sphere vertex expressed as a (lon, lat) tuple given in degrees
    pub fn proj_lonlat(&self, lonlat: &LonLat) -> Option<ImgXY> {
        let img_xy = match self {
            // Zenithal
            WCSProj::Azp(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Szp(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Tan(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Stg(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Sin(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Arc(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Zpn(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Zea(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Air(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Ncp(wcs) => wcs.lonlat2img(lonlat),
            // Pseudo-cyl
            WCSProj::Cyp(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Cea(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Car(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Mer(wcs) => wcs.lonlat2img(lonlat),
            // Cylindrical
            WCSProj::Sfl(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Par(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Mol(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Ait(wcs) => wcs.lonlat2img(lonlat),
            // Conic
            WCSProj::Cop(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Cod(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Coe(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::Coo(wcs) => wcs.lonlat2img(lonlat),

            /* Sip variants */
            // Zenithal
            WCSProj::AzpSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::SzpSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::TanSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::StgSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::SinSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::ArcSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::ZpnSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::ZeaSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::AirSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::NcpSip(wcs) => wcs.lonlat2img(lonlat),
            // Pseudo-cyl
            WCSProj::CypSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::CeaSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::CarSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::MerSip(wcs) => wcs.lonlat2img(lonlat),
            // Cylindrical
            WCSProj::SflSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::ParSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::MolSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::AitSip(wcs) => wcs.lonlat2img(lonlat),
            // Conic
            WCSProj::CopSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::CodSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::CoeSip(wcs) => wcs.lonlat2img(lonlat),
            WCSProj::CooSip(wcs) => wcs.lonlat2img(lonlat),
        };

        img_xy.map(|xy| ImgXY::new(xy.x() - 1.0, xy.y() - 1.0))
    }

    /// Unproject a (X, Y) point from the image space to get its corresponding location on the sphere
    /// The result is given a (lon, lat) tuple expressed in degrees.
    /// 
    /// # Param
    /// * `img_pos`: the image space point expressed as a (X, Y) tuple given en pixels
    pub fn unproj_lonlat(&self, img_pos: &ImgXY) -> Option<LonLat> {
        let img_pos = ImgXY::new(img_pos.x() + 1.0, img_pos.y() + 1.0);
        match self {
            // Zenithal
            WCSProj::Azp(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Szp(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Tan(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Stg(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Sin(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Arc(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Zpn(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Zea(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Air(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Ncp(wcs) => wcs.img2lonlat(&img_pos),
            // Pseudo-cyl
            WCSProj::Cyp(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Cea(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Car(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Mer(wcs) => wcs.img2lonlat(&img_pos),
            // Cylindrical
            WCSProj::Sfl(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Par(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Mol(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Ait(wcs) => wcs.img2lonlat(&img_pos),
            // Conic
            WCSProj::Cop(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Cod(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Coe(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::Coo(wcs) => wcs.img2lonlat(&img_pos),

            /* Sip variants */
            // Zenithal
            WCSProj::AzpSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::SzpSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::TanSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::StgSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::SinSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::ArcSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::ZpnSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::ZeaSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::AirSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::NcpSip(wcs) => wcs.img2lonlat(&img_pos),
            // Pseudo-cyl
            WCSProj::CypSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::CeaSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::CarSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::MerSip(wcs) => wcs.img2lonlat(&img_pos),
            // Cylindrical
            WCSProj::SflSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::ParSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::MolSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::AitSip(wcs) => wcs.img2lonlat(&img_pos),
            // Conic
            WCSProj::CopSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::CodSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::CoeSip(wcs) => wcs.img2lonlat(&img_pos),
            WCSProj::CooSip(wcs) => wcs.img2lonlat(&img_pos),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::WCS;
    use crate::mapproj::Projection;
    use fitsrs::fits::Fits;
    use fitsrs::hdu::{
        header::{
            extension::image::Image,
            Header,
        },
        data::iter
    };
    use glob::glob;
    use mapproj::{CanonicalProjection, ImgXY, LonLat};
    use std::f64::consts::PI;
    use std::fs::File;
    use std::io::BufReader;

    #[test]
    fn test_visualize() {
        let f = File::open("examples/panstarrs-rotated-around-orion.fits").unwrap();

        let mut reader = BufReader::new(f);
        let Fits { mut hdu} = Fits::from_reader(&mut reader).unwrap();

        // Parse data
        let data = match hdu.get_data_mut() {
            iter::Data::F32(it) => it.collect::<Vec<_>>(),
            _ => unreachable!(),
        };

        let header = hdu.get_header();
        let wcs = WCS::new(&header).unwrap();
        reproject_fits_image(mapproj::zenithal::azp::Azp::new(), &wcs, &header, &data);
        reproject_fits_image(mapproj::zenithal::szp::Szp::new(), &wcs, &header, &data);
        reproject_fits_image(mapproj::zenithal::tan::Tan::new(), &wcs, &header, &data);
        reproject_fits_image(mapproj::zenithal::stg::Stg::new(), &wcs, &header, &data);
        reproject_fits_image(mapproj::zenithal::sin::Sin::new(), &wcs, &header, &data);
        reproject_fits_image(mapproj::zenithal::arc::Arc::new(), &wcs, &header, &data);
        reproject_fits_image(mapproj::zenithal::zea::Zea::new(), &wcs, &header, &data);
        reproject_fits_image(mapproj::zenithal::air::Air::new(), &wcs, &header, &data);

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
    }

    fn reproject_fits_image<'a, T: CanonicalProjection>(
        proj: T,
        wcs: &WCS,
        header: &Header<Image>,
        data: &[f32],
    ) {
        let scale = header
            .get_parsed::<f64>(b"BSCALE  ")
            .unwrap_or(Ok(1.0))
            .unwrap() as f32;
        let offset = header
            .get_parsed::<f64>(b"BZERO   ")
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
        //let f = File::open("examples/cutout-CDS_P_HST_PHAT_F475W.fits").unwrap();
        let f = File::open("examples/cutout-CDS_P_PanSTARRS_DR1_g.fits").unwrap();
        //let f = File::open("examples/FOCx38i0101t_c0f.fits").unwrap();
        //let f = File::open("examples/pc.fits").unwrap();
        let mut reader = BufReader::new(f);
        let Fits { hdu } = Fits::from_reader(&mut reader).unwrap();
        let header = hdu.get_header();
        let wcs = WCS::new(header).unwrap();

        use std::fs::File;
        // Build the CSV reader and iterate over each record.
        //let f = File::open("examples/pc.fits.csv").unwrap();
        let f = File::open("examples/cutout-CDS_P_PanSTARRS_DR1_g.fits.csv").unwrap();
        //let f = File::open("examples/FOCx38i0101t_c0f.fits.csv").unwrap();
        //let f = File::open("examples/cutout-CDS_P_HST_PHAT_F475W.fits.csv").unwrap();

        let mut rdr = csv::Reader::from_reader(BufReader::new(f));
        for result in rdr.records() {
            let record = result.unwrap();

            let ra: f64 = record[0].parse().unwrap();
            let dec: f64 = record[1].parse().unwrap();
            let x: f64 = record[2].parse().unwrap();
            let y: f64 = record[3].parse().unwrap();

            if ra.is_finite() && dec.is_finite() {
                if let Some(img_xy) = wcs.proj(&LonLat::new(ra, dec)) {
                    assert_delta!(img_xy.x(), x, 1e-4);
                    assert_delta!(img_xy.y(), y, 1e-4);
                }
            }
        }
    }

    #[test]
    fn crval_to_crpix() {
        for entry in glob("examples/*.fits").unwrap() {
            if let Ok(path) = dbg!(entry) {
                let f = File::open(path).unwrap();
                let mut reader = BufReader::new(f);
                let Fits { hdu} = Fits::from_reader(&mut reader).unwrap();
                let header = hdu.get_header();
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

                let wcs = WCS::new(&header).unwrap();

                // crval to crpix
                let proj_px = wcs
                    .proj(&LonLat::new(dbg!(crval1).to_radians(), dbg!(crval2).to_radians()))
                    .unwrap();
                assert_delta!(proj_px.x(), crpix1 - 1.0, 1e-6);
                assert_delta!(proj_px.y(), crpix2 - 1.0, 1e-6);

                // crpix to crval
                let lonlat = wcs
                    .unproj_lonlat(&ImgXY::new(dbg!(crpix1) - 1.0, dbg!(crpix2) - 1.0))
                    .unwrap();
                assert_delta!(lonlat.lon(), crval1.to_radians(), 1e-6);
                assert_delta!(lonlat.lat(), crval2.to_radians(), 1e-6);
            }
        }
    }

    /*#[test]
    fn open_fits() {
        let f = File::open("examples/M81_B_GaiaHduImg1.fits").unwrap();

        let mut reader = BufReader::new(f);
        let Fits { hdu} = Fits::from_reader(&mut reader).unwrap();
        let header = hdu.get_header();
        assert!(WCS::new(header).is_err());
    }*/
}
