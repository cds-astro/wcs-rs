extern crate fitsrs;
extern crate mapproj;
#[macro_use]
extern crate quick_error;

mod error;
use error::Error;
mod projection;

use crate::projection::WCSCanonicalProjection;

// Imports
use fitsrs::hdu::{DataRead, HDU};
use mapproj::{
    img2celestial::Img2Celestial,
    img2proj::WcsImgXY2ProjXY,
    LonLat,
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
}

pub type ImgXY = mapproj::ImgXY;
pub type XYZ = mapproj::XYZ;

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
        let ctype1 = header
            .get_parsed::<String>(b"CTYPE1  ")
            .unwrap_or(Ok(String::from("RA---TAN")))?;
        let _ = header
            .get_parsed::<String>(b"CTYPE2  ")
            .unwrap_or(Ok(String::from("DEC--TAN")))?;

        let proj_name = &ctype1[5..];
        match proj_name {
            // Zenithal
            "AZP" => {
                let proj = Azp::parse_proj(&header)?;
                Ok(WCS::Azp(Img2Celestial::new(img2proj, proj)))
            }
            "SZP" => {
                let proj = Szp::parse_proj(&header)?;
                Ok(WCS::Szp(Img2Celestial::new(img2proj, proj)))
            }
            "TAN" => {
                let proj = Tan::parse_proj(&header)?;
                Ok(WCS::Tan(Img2Celestial::new(img2proj, proj)))
            }
            "STG" => {
                let proj = Stg::parse_proj(&header)?;
                Ok(WCS::Stg(Img2Celestial::new(img2proj, proj)))
            }
            "SIN" => {
                let proj = Sin::parse_proj(&header)?;
                Ok(WCS::Sin(Img2Celestial::new(img2proj, proj)))
            }
            "ARC" => {
                let proj = Arc::parse_proj(&header)?;
                Ok(WCS::Arc(Img2Celestial::new(img2proj, proj)))
            }
            "ZPN" => {
                let proj = Zpn::parse_proj(&header)?;
                Ok(WCS::Zpn(Img2Celestial::new(img2proj, proj)))
            }
            "ZEA" => {
                let proj = Zea::parse_proj(&header)?;
                Ok(WCS::Zea(Img2Celestial::new(img2proj, proj)))
            }
            "AIR" => {
                let proj = Air::parse_proj(&header)?;
                Ok(WCS::Air(Img2Celestial::new(img2proj, proj)))
            }
            // Cylindrical
            "CYP" => {
                let proj = Cyp::parse_proj(&header)?;
                Ok(WCS::Cyp(Img2Celestial::new(img2proj, proj)))
            }
            "CEA" => {
                let proj = Cea::parse_proj(&header)?;
                Ok(WCS::Cea(Img2Celestial::new(img2proj, proj)))
            }
            "CAR" => {
                let proj = Car::parse_proj(&header)?;
                Ok(WCS::Car(Img2Celestial::new(img2proj, proj)))
            }
            "MER" => {
                let proj = Mer::parse_proj(&header)?;
                Ok(WCS::Mer(Img2Celestial::new(img2proj, proj)))
            }
            // Pseudo-cylindrical
            "SFL" => {
                let proj = Sfl::parse_proj(&header)?;
                Ok(WCS::Sfl(Img2Celestial::new(img2proj, proj)))
            }
            "PAR" => {
                let proj = Par::parse_proj(&header)?;
                Ok(WCS::Par(Img2Celestial::new(img2proj, proj)))
            }
            "MOL" => {
                let proj = Mol::parse_proj(&header)?;
                Ok(WCS::Mol(Img2Celestial::new(img2proj, proj)))
            }
            "AIT" => {
                let proj = Ait::parse_proj(&header)?;
                Ok(WCS::Ait(Img2Celestial::new(img2proj, proj)))
            }
            // Conic
            "COP" => {
                let proj = Cop::parse_proj(&header)?;
                Ok(WCS::Cop(Img2Celestial::new(img2proj, proj)))
            }
            "COD" => {
                let proj = Cod::parse_proj(&header)?;
                Ok(WCS::Cod(Img2Celestial::new(img2proj, proj)))
            }
            "COE" => {
                let proj = Coe::parse_proj(&header)?;
                Ok(WCS::Coe(Img2Celestial::new(img2proj, proj)))
            }
            "COO" => {
                let proj = Coo::parse_proj(&header)?;
                Ok(WCS::Coo(Img2Celestial::new(img2proj, proj)))
            }
            _ => Err(Error::NotImplementedProjection(proj_name.to_string())),
        }
    }

    /// Project a 3D vector (XYZ) to the image space (xy) in pixels
    ///
    ///
    pub fn proj(&self, xyz: &XYZ) -> Option<ImgXY> {
        let img_xy = match self {
            // Zenithal
            WCS::Azp(wcs) => wcs.xyz2img(xyz),
            WCS::Szp(wcs) => wcs.xyz2img(xyz),
            WCS::Tan(wcs) => wcs.xyz2img(xyz),
            WCS::Stg(wcs) => wcs.xyz2img(xyz),
            WCS::Sin(wcs) => wcs.xyz2img(xyz),
            WCS::Arc(wcs) => wcs.xyz2img(xyz),
            WCS::Zpn(wcs) => wcs.xyz2img(xyz),
            WCS::Zea(wcs) => wcs.xyz2img(xyz),
            WCS::Air(wcs) => wcs.xyz2img(xyz),
            // Pseudo-cyl
            WCS::Cyp(wcs) => wcs.xyz2img(xyz),
            WCS::Cea(wcs) => wcs.xyz2img(xyz),
            WCS::Car(wcs) => wcs.xyz2img(xyz),
            WCS::Mer(wcs) => wcs.xyz2img(xyz),
            // Cylindrical
            WCS::Sfl(wcs) => wcs.xyz2img(xyz),
            WCS::Par(wcs) => wcs.xyz2img(xyz),
            WCS::Mol(wcs) => wcs.xyz2img(xyz),
            WCS::Ait(wcs) => wcs.xyz2img(xyz),
            // Conic
            WCS::Cop(wcs) => wcs.xyz2img(xyz),
            WCS::Cod(wcs) => wcs.xyz2img(xyz),
            WCS::Coe(wcs) => wcs.xyz2img(xyz),
            WCS::Coo(wcs) => wcs.xyz2img(xyz),
        };

        img_xy.map(|xy| ImgXY::new(xy.x() - 1.0, xy.y() - 1.0))
    }

    pub fn proj_lonlat(&self, lon: f64, lat: f64) -> Option<ImgXY> {
        let xyz = &LonLat::new(lon, lat).to_xyz();
        let img_xy = match self {
            // Zenithal
            WCS::Azp(wcs) => wcs.xyz2img(xyz),
            WCS::Szp(wcs) => wcs.xyz2img(xyz),
            WCS::Tan(wcs) => wcs.xyz2img(xyz),
            WCS::Stg(wcs) => wcs.xyz2img(xyz),
            WCS::Sin(wcs) => wcs.xyz2img(xyz),
            WCS::Arc(wcs) => wcs.xyz2img(xyz),
            WCS::Zpn(wcs) => wcs.xyz2img(xyz),
            WCS::Zea(wcs) => wcs.xyz2img(xyz),
            WCS::Air(wcs) => wcs.xyz2img(xyz),
            // Pseudo-cyl
            WCS::Cyp(wcs) => wcs.xyz2img(xyz),
            WCS::Cea(wcs) => wcs.xyz2img(xyz),
            WCS::Car(wcs) => wcs.xyz2img(xyz),
            WCS::Mer(wcs) => wcs.xyz2img(xyz),
            // Cylindrical
            WCS::Sfl(wcs) => wcs.xyz2img(xyz),
            WCS::Par(wcs) => wcs.xyz2img(xyz),
            WCS::Mol(wcs) => wcs.xyz2img(xyz),
            WCS::Ait(wcs) => wcs.xyz2img(xyz),
            // Conic
            WCS::Cop(wcs) => wcs.xyz2img(xyz),
            WCS::Cod(wcs) => wcs.xyz2img(xyz),
            WCS::Coe(wcs) => wcs.xyz2img(xyz),
            WCS::Coo(wcs) => wcs.xyz2img(xyz),
        };

        img_xy.map(|xy| ImgXY::new(xy.x() - 1.0, xy.y() - 1.0))
    }

    pub fn unproj(&self, img_pos: &ImgXY) -> Option<XYZ> {
        let img_pos = ImgXY::new(img_pos.x() + 1.0, img_pos.y() + 1.0);
        match self {
            // Zenithal
            WCS::Azp(wcs) => wcs.img2xyz(&img_pos),
            WCS::Szp(wcs) => wcs.img2xyz(&img_pos),
            WCS::Tan(wcs) => wcs.img2xyz(&img_pos),
            WCS::Stg(wcs) => wcs.img2xyz(&img_pos),
            WCS::Sin(wcs) => wcs.img2xyz(&img_pos),
            WCS::Arc(wcs) => wcs.img2xyz(&img_pos),
            WCS::Zpn(wcs) => wcs.img2xyz(&img_pos),
            WCS::Zea(wcs) => wcs.img2xyz(&img_pos),
            WCS::Air(wcs) => wcs.img2xyz(&img_pos),
            // Pseudo-cyl
            WCS::Cyp(wcs) => wcs.img2xyz(&img_pos),
            WCS::Cea(wcs) => wcs.img2xyz(&img_pos),
            WCS::Car(wcs) => wcs.img2xyz(&img_pos),
            WCS::Mer(wcs) => wcs.img2xyz(&img_pos),
            // Cylindrical
            WCS::Sfl(wcs) => wcs.img2xyz(&img_pos),
            WCS::Par(wcs) => wcs.img2xyz(&img_pos),
            WCS::Mol(wcs) => wcs.img2xyz(&img_pos),
            WCS::Ait(wcs) => wcs.img2xyz(&img_pos),
            // Conic
            WCS::Cop(wcs) => wcs.img2xyz(&img_pos),
            WCS::Cod(wcs) => wcs.img2xyz(&img_pos),
            WCS::Coe(wcs) => wcs.img2xyz(&img_pos),
            WCS::Coo(wcs) => wcs.img2xyz(&img_pos),
        }
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
        }
    }
}

#[cfg(test)]
mod tests {
    use super::WCS;
    use fitsrs::fits::Fits;
    use fitsrs::hdu::data::DataOwned;
    use mapproj::CanonicalProjection;
    use mapproj::{ImgXY, LonLat};
    use std::fs::File;
    use std::io::BufReader;

    use glob::glob;

    #[test]
    fn reproject_fits_image() {
        //let f = File::open("examples/FOCx38i0101t_c0f.fits").unwrap();
        let f = File::open("examples/cutout-CDS_P_PanSTARRS_DR1_g.fits").unwrap();

        let Fits { hdu } = Fits::from_reader(BufReader::new(f)).unwrap();

        let wcs = WCS::new(&hdu).unwrap();

        let header = &hdu.header;
        let scale = header
            .get_parsed::<f64>(b"BSCALE  ")
            .unwrap_or(Ok(1.0))
            .unwrap() as f32;
        let offset = header
            .get_parsed::<f64>(b"BZERO   ")
            .unwrap_or(Ok(0.0))
            .unwrap() as f32;

        let data = hdu.data;

        let width = *header.get_axis_size(1).unwrap();
        let height = *header.get_axis_size(2).unwrap();

        // Parse data
        let data = match data {
            DataOwned::F32(it) => it.collect::<Vec<_>>(),
            _ => unreachable!(),
        };

        // Define mollweide projection in which we will project
        let proj = mapproj::pseudocyl::mol::Mol::new();
        //let proj = mapproj::zenithal::sin::Sin::new();
        let bounds = proj.bounds();
        let x_bounds = bounds.x_bounds().as_ref().unwrap();
        let y_bounds = bounds.y_bounds().as_ref().unwrap();

        let x_off = x_bounds.start();
        let x_len = x_bounds.end() - x_bounds.start();

        let y_off = y_bounds.start();
        let y_len = y_bounds.end() - y_bounds.start();

        const WIDTH_IMAGE: usize = 1024;
        const HEIGHT_IMAGE: usize = WIDTH_IMAGE / 2;

        // Create a new ImgBuf with width: imgx and height: imgy
        let mut imgbuf = image::ImageBuffer::new(WIDTH_IMAGE as u32, HEIGHT_IMAGE as u32);

        for y in 0..height {
            for x in 0..width {
                let grayscale_val = (data[y * width + x] * scale + offset) as u8;

                let img_xy = ImgXY::new(x as f64, y as f64);
                if let Some(xyz) = wcs.unproj(&img_xy) {
                    if let Some(proj_xy) = proj.proj(&xyz) {
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

        // Save the image as “fractal.png”, the format is deduced from the path
        imgbuf.save("fits-viewer.jpeg").unwrap();
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
        //let f = File::open("examples/cutout-CDS_P_PanSTARRS_DR1_g.fits").unwrap();
        //let f = File::open("examples/FOCx38i0101t_c0f.fits").unwrap();
        let f = File::open("examples/pc.fits").unwrap();

        let Fits { hdu } = Fits::from_reader(BufReader::new(f)).unwrap();
        let wcs = WCS::new(&hdu).unwrap();

        use std::fs::File;
        // Build the CSV reader and iterate over each record.
        let f = File::open("examples/pc.fits.csv").unwrap();
        //let f = File::open("examples/cutout-CDS_P_PanSTARRS_DR1_g.fits.csv").unwrap();
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
                if let Some(img_xy) = wcs.proj_lonlat(ra, dec) {
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
                    .proj_lonlat(crval1.to_radians(), crval2.to_radians())
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
