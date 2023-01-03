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
use fitsrs::hdu::{DataRead, HDU};
use mapproj::{
    img2celestial::Img2Celestial,
    img2proj::{WcsImgXY2ProjXY, WcsWithSipImgXY2ProjXY},
    conic::{cod::Cod, coe::Coe, coo::Coo, cop::Cop},
    cylindrical::{car::Car, cea::Cea, cyp::Cyp, mer::Mer},
    pseudocyl::{ait::Ait, mol::Mol, par::Par, sfl::Sfl},
    zenithal::{
        air::Air,
        arc::Arc,
        azp::Azp,
        sin::Sin,
        stg::Stg,
        szp::Szp,
        tan::Tan,
        zea::Zea,
        zpn::Zpn,
    },
};
use paste::paste;
macro_rules! create_specific_proj {
    ( $proj_name:ident, $header:expr, $ctype1:expr, $crpix1:expr, $crpix2:expr, $img2proj:expr ) => {
        {
            let proj = $proj_name::parse_proj(&$header)?;

            let is_sip_found = &$ctype1[($ctype1.len() - 3)..] == "SIP";
            if is_sip_found {
                let sip = sip::parse_sip($header, $crpix1, $crpix2)?;
                let img2proj = WcsWithSipImgXY2ProjXY::new($img2proj, sip);

                paste! {    
                    Ok(WCS::[ <$proj_name Sip> ](Img2Celestial::new(img2proj, proj)))
                }
            } else {
                Ok(WCS::$proj_name(Img2Celestial::new($img2proj, proj)))
            }
        }
    };
}


pub enum WCS {
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

pub type ImgXY = mapproj::ImgXY;
pub type LonLat = mapproj::LonLat;

impl WCS {
    pub fn new<'a, R>(hdu: &HDU<'a, R>) -> Result<Self, Error>
    where
        R: DataRead<'a> + 'a,
    {
        let header = &hdu.header;
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

            WcsImgXY2ProjXY::from_cd(
                crpix1,
                crpix2,
                cd11,
                cd12,
                cd21,
                cd22,
            )
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

    pub fn proj_lonlat(&self, lonlat: &LonLat) -> Option<ImgXY> {
        let img_xy = match self {
            // Zenithal
            WCS::Azp(wcs) => wcs.lonlat2img(lonlat),
            WCS::Szp(wcs) => wcs.lonlat2img(lonlat),
            WCS::Tan(wcs) => wcs.lonlat2img(lonlat),
            WCS::Stg(wcs) => wcs.lonlat2img(lonlat),
            WCS::Sin(wcs) => wcs.lonlat2img(lonlat),
            WCS::Arc(wcs) => wcs.lonlat2img(lonlat),
            WCS::Zpn(wcs) => wcs.lonlat2img(lonlat),
            WCS::Zea(wcs) => wcs.lonlat2img(lonlat),
            WCS::Air(wcs) => wcs.lonlat2img(lonlat),
            // Pseudo-cyl
            WCS::Cyp(wcs) => wcs.lonlat2img(lonlat),
            WCS::Cea(wcs) => wcs.lonlat2img(lonlat),
            WCS::Car(wcs) => wcs.lonlat2img(lonlat),
            WCS::Mer(wcs) => wcs.lonlat2img(lonlat),
            // Cylindrical
            WCS::Sfl(wcs) => wcs.lonlat2img(lonlat),
            WCS::Par(wcs) => wcs.lonlat2img(lonlat),
            WCS::Mol(wcs) => wcs.lonlat2img(lonlat),
            WCS::Ait(wcs) => wcs.lonlat2img(lonlat),
            // Conic
            WCS::Cop(wcs) => wcs.lonlat2img(lonlat),
            WCS::Cod(wcs) => wcs.lonlat2img(lonlat),
            WCS::Coe(wcs) => wcs.lonlat2img(lonlat),
            WCS::Coo(wcs) => wcs.lonlat2img(lonlat),

            /* Sip variants */

            // Zenithal
            WCS::AzpSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::SzpSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::TanSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::StgSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::SinSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::ArcSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::ZpnSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::ZeaSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::AirSip(wcs) => wcs.lonlat2img(lonlat),
            // Pseudo-cyl
            WCS::CypSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::CeaSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::CarSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::MerSip(wcs) => wcs.lonlat2img(lonlat),
            // Cylindrical
            WCS::SflSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::ParSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::MolSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::AitSip(wcs) => wcs.lonlat2img(lonlat),
            // Conic
            WCS::CopSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::CodSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::CoeSip(wcs) => wcs.lonlat2img(lonlat),
            WCS::CooSip(wcs) => wcs.lonlat2img(lonlat),
        };

        img_xy.map(|xy| ImgXY::new(xy.x() - 1.0, xy.y() - 1.0))
    }

    pub fn unproj_lonlat(&self, img_pos: &ImgXY) -> Option<LonLat> {
        let img_pos = ImgXY::new(img_pos.x() + 1.0, img_pos.y() + 1.0);
        match self {
            // Zenithal
            WCS::Azp(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Szp(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Tan(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Stg(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Sin(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Arc(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Zpn(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Zea(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Air(wcs) => wcs.img2lonlat(&img_pos),
            // Pseudo-cyl
            WCS::Cyp(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Cea(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Car(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Mer(wcs) => wcs.img2lonlat(&img_pos),
            // Cylindrical
            WCS::Sfl(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Par(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Mol(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Ait(wcs) => wcs.img2lonlat(&img_pos),
            // Conic
            WCS::Cop(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Cod(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Coe(wcs) => wcs.img2lonlat(&img_pos),
            WCS::Coo(wcs) => wcs.img2lonlat(&img_pos),

            /* Sip variants */
            // Zenithal
            WCS::AzpSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::SzpSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::TanSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::StgSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::SinSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::ArcSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::ZpnSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::ZeaSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::AirSip(wcs) => wcs.img2lonlat(&img_pos),
            // Pseudo-cyl
            WCS::CypSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::CeaSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::CarSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::MerSip(wcs) => wcs.img2lonlat(&img_pos),
            // Cylindrical
            WCS::SflSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::ParSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::MolSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::AitSip(wcs) => wcs.img2lonlat(&img_pos),
            // Conic
            WCS::CopSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::CodSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::CoeSip(wcs) => wcs.img2lonlat(&img_pos),
            WCS::CooSip(wcs) => wcs.img2lonlat(&img_pos),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::WCS;
    use fitsrs::fits::Fits;
    use fitsrs::hdu::{HDU, Header};
    use fitsrs::hdu::data::DataOwned;
    use mapproj::{ImgXY, LonLat, CanonicalProjection};
    use std::f64::consts::PI;
    use std::fs::File;
    use std::io::BufReader;
    use crate::mapproj::Projection;
    use glob::glob;

    #[test]
    fn test_visualize() {
        let f = File::open("examples/cutout-CDS_P_PanSTARRS_DR1_g.fits").unwrap();

        let Fits { hdu } = Fits::from_reader(BufReader::new(f)).unwrap();
        let wcs = WCS::new(&hdu).unwrap();

        // Parse data
        let HDU { header, data } = hdu;
        let data = match data {
            DataOwned::F32(it) => it.collect::<Vec<_>>(),
            _ => unreachable!(),
        };

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

    fn reproject_fits_image<'a, T: CanonicalProjection>(proj: T, wcs: &WCS, header: &Header, data: &[f32]) {
        let scale = header
            .get_parsed::<f64>(b"BSCALE  ")
            .unwrap_or(Ok(1.0))
            .unwrap() as f32;
        let offset = header
            .get_parsed::<f64>(b"BZERO   ")
            .unwrap_or(Ok(0.0))
            .unwrap() as f32;


        let width = *header.get_axis_size(1).unwrap();
        let height = *header.get_axis_size(2).unwrap();

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
                let grayscale_val = (data[y * width + x] * scale + offset) as u8;

                let img_xy = ImgXY::new(x as f64, y as f64);
                if let Some(lonlat) = wcs.unproj_lonlat(&img_xy) {
                    if let Some(proj_xy) = proj.proj_lonlat(&lonlat) {
                        let proj_x = ((proj_xy.x() as f64) - x_off) / x_len; // between 0 and 1
                        let proj_y = ((proj_xy.y() as f64) - y_off) / y_len; // between 0 and 1

                        if (0.0..1.0).contains(&proj_x) && (0.0..1.0).contains(&proj_y) {
                            let ix = (proj_x * (WIDTH_IMAGE as f64)) as usize;
                            let iy = (proj_y * (HEIGHT_IMAGE as f64)) as usize;

                            let pixel = imgbuf.get_pixel_mut(ix as u32, iy as u32);
                            *pixel = image::Rgb([grayscale_val, grayscale_val, grayscale_val]);
                        }
                    }
                }
            }
        }

        let filename = &format!("tests/reproj/fits-{}.jpeg", <T as CanonicalProjection>::WCS_NAME);
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

        let Fits { hdu } = Fits::from_reader(BufReader::new(f)).unwrap();
        let wcs = WCS::new(&hdu).unwrap();

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
                if let Some(img_xy) = wcs.proj_lonlat(&LonLat::new(ra, dec)) {
                    println!("{:?}", img_xy.x() - x);
                    println!("{:?}", img_xy.y() - y);

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

                let Fits { hdu } = Fits::from_reader(BufReader::new(f)).unwrap();
                let header = &hdu.header;
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
                let wcs = WCS::new(&hdu).unwrap();

                // crval to crpix
                let proj_px = wcs
                    .proj_lonlat(&LonLat::new(crval1.to_radians(), crval2.to_radians()))
                    .unwrap();
                assert_delta!(proj_px.x(), crpix1 - 1.0, 1e-6);
                assert_delta!(proj_px.y(), crpix2 - 1.0, 1e-6);

                // crpix to crval
                let lonlat = wcs
                    .unproj_lonlat(&ImgXY::new(crpix1 - 1.0, crpix2 - 1.0))
                    .unwrap();
                assert_delta!(lonlat.lon(), crval1.to_radians(), 1e-6);
                assert_delta!(lonlat.lat(), crval2.to_radians(), 1e-6);
            }
        }
    }
}
