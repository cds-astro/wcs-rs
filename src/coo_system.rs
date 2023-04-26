use fitsrs::hdu::header::{
    Header,
    extension::image::Image
};

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
    GAPPT
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
            _ => {
                Err(Error::UnrecognizedRadeSys(radesys))
            }
        }
    }
}


pub enum CooSystem {
    EQUATORIAL,
    GALACTIC,
    ECLIPTIC,
    HELIOECLIPTIC,
    SUPERGALACTIC,
    CUSTOM {
        radesys: RadeSys,
        equinox: f64,
    }
}

impl CooSystem {
    pub fn parse(header: &Header<Image>) -> Result<Self, Error> {
        let ctype1 = utils::retrieve_mandatory_parsed_keyword::<String>(header, "CTYPE1  ")?;

        let coo_system = match &ctype1[0..4] {
            "RA--" => CooSystem::EQUATORIAL,
            "GLON" => CooSystem::GALACTIC,
            "ELON" => CooSystem::ECLIPTIC,
            "HLON" => CooSystem::HELIOECLIPTIC,
            "SLON" => CooSystem::SUPERGALACTIC,
            _ => {
                let radesys = RadeSys::parse(header)?;
                let equinox = utils::retrieve_mandatory_parsed_keyword::<f64>(header, "EQUINOX ")?;

                CooSystem::CUSTOM { radesys, equinox }
            }
        };

        Ok(coo_system)
    }
}