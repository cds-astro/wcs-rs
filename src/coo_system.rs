use std::f64::consts::PI;

use fitsrs::hdu::header::{extension::image::Image, Header};
use mapproj::LonLat;
use mapproj::XYZ;

use crate::error::Error;
use crate::utils;

pub enum RadeSys {
    /// International Celestial Reference System
    ICRS,
    /// mean place, new (IAU 1984) system
    Fk5,
    /// mean place, old (Bessell-Newcomb) system
    Fk4,
    /// mean place, old system but without e-terms
    Fk4NoE,
    /// geocentric apparent place, IAU 1984 system
    GAPPT,
}

impl RadeSys {
    pub fn parse(header: &Header<Image>) -> Result<Self, Error> {
        let radesys = utils::retrieve_mandatory_parsed_keyword::<String>(header, "RADESYS ")?;

        match radesys.as_str() {
            "ICRS" => Ok(RadeSys::ICRS),
            "FK5" => Ok(RadeSys::Fk5),
            "FK4" => Ok(RadeSys::Fk4),
            "FK4-NO-E" => Ok(RadeSys::Fk4NoE),
            "GAPPT" => Ok(RadeSys::GAPPT),
            _ => Err(Error::UnrecognizedRadeSys(radesys)),
        }
    }
}

pub enum CooSystem {
    /// ICRS/J2000
    EQUATORIAL,
    GALACTIC,
    ECLIPTIC,
    HELIOECLIPTIC,
    SUPERGALACTIC,
    CUSTOM {
        radesys: RadeSys,
        equinox: f64,
    },
}

impl CooSystem {
    pub fn parse(header: &Header<Image>) -> Result<Self, Error> {
        let equinox = utils::retrieve_mandatory_parsed_keyword::<f64>(header, "EQUINOX ");
        let radesys = RadeSys::parse(header);

        let coo_system = if let (Ok(radesys), Ok(equinox)) = (radesys, equinox) {
            // if there is a radesys take it into account
            CooSystem::CUSTOM { radesys, equinox }
        } else {
            let ctype1 = utils::retrieve_mandatory_parsed_keyword::<String>(header, "CTYPE1  ")?;

            match ctype1.as_bytes()[0] {
                b'G' => CooSystem::GALACTIC,
                b'E' => CooSystem::ECLIPTIC,
                b'H' => CooSystem::HELIOECLIPTIC,
                b'S' => CooSystem::SUPERGALACTIC,
                _ => CooSystem::EQUATORIAL,
            }
        };

        Ok(coo_system)
    }

    /// Convert a lonlat expressed in icrs coordinate system to the
    /// self system
    pub fn from_icrs(&self, lonlat: LonLat) -> LonLat {
        match self {
            CooSystem::EQUATORIAL => lonlat,
            CooSystem::GALACTIC => {
                let xyz = lonlat.to_xyz();
                // ICRS_2_GAL * xyz
                let rotated_xyz = XYZ::new(
                    ICRS_2_GAL[0] * xyz.x() + ICRS_2_GAL[1] * xyz.y() + ICRS_2_GAL[2] * xyz.z(),
                    ICRS_2_GAL[3] * xyz.x() + ICRS_2_GAL[4] * xyz.y() + ICRS_2_GAL[5] * xyz.z(),
                    ICRS_2_GAL[6] * xyz.x() + ICRS_2_GAL[7] * xyz.y() + ICRS_2_GAL[8] * xyz.z(),
                );

                xyz_to_lonlat(&rotated_xyz)
            }
            _ => {
                // todo
                lonlat
            }
        }
    }

    /// Convert a lonlat expressed in the self system to the icrs coo system
    pub fn to_icrs(&self, lonlat: LonLat) -> LonLat {
        match self {
            CooSystem::EQUATORIAL => lonlat,
            CooSystem::GALACTIC => {
                let xyz = lonlat.to_xyz();
                // ICRS_2_GAL * xyz
                let rotated_xyz = XYZ::new(
                    GAL_2_ICRS[0] * xyz.x() + GAL_2_ICRS[1] * xyz.y() + GAL_2_ICRS[2] * xyz.z(),
                    GAL_2_ICRS[3] * xyz.x() + GAL_2_ICRS[4] * xyz.y() + GAL_2_ICRS[5] * xyz.z(),
                    GAL_2_ICRS[6] * xyz.x() + GAL_2_ICRS[7] * xyz.y() + GAL_2_ICRS[8] * xyz.z(),
                );

                xyz_to_lonlat(&rotated_xyz)
            }
            _ => {
                // todo
                lonlat
            }
        }
    }
}

const GAL_2_ICRS: &[f64; 9] = &[
    -0.4448296299195045,
    -0.1980763734646737,
    -0.873437090247923,
    0.7469822444763707,
    0.4559837762325372,
    -0.4838350155267381,
    0.4941094279435681,
    -0.867_666_148_981_161,
    -0.0548755604024359,
];

const ICRS_2_GAL: &[f64; 9] = &[
    -0.4448296299195045,
    0.7469822444763707,
    0.4941094279435681,
    -0.1980763734646737,
    0.4559837762325372,
    -0.867_666_148_981_161,
    -0.873437090247923,
    -0.4838350155267381,
    -0.0548755604024359,
];

use mapproj::CustomFloat;
fn xyz_to_lonlat(xyz: &XYZ) -> LonLat {
    let r2 = xyz.x().pow2() + xyz.y().pow2();
    // Latitude in [-pi/2, pi/2] (ok, since cos always positive here)
    let lat = xyz.z().atan2(r2.sqrt());
    // Compute the longitude in [-pi, pi]
    let lon = xyz.y().atan2(xyz.x());
    // Conforms to convention: Longitude in [0, 2*PI]
    LonLat::new(if lon < 0.0 { 2.0 * PI + lon } else { lon }, lat)
}
