quick_error! {
    #[derive(Debug)]
    pub enum Error {
        IoFile(err: std::io::Error) {
            from()
        }
        /// Encountered when parsing the Fits file, returned from fitsrs
        CardWrongType(key: String, tt: String) {
            display("{} could not be parsed {}", key, tt)
        }
        /// Encountered when parsing the Fits file, returned from fitsrs
        FitsParsing(err: fitsrs::error::Error) {
            from()
            display("{}", err)
        }
        /// Unprojection not defined
        UnprojNotDefined(img_x: f64, img_y: f64) {
            display("({}, {}) is out of projection", img_x, img_y)
        }
        /// Encountered when the projection cannot be instantiated
        /// i.e. Zpn needs coefficients meeting a specific condition
        InitProjection(proj_name: &'static str, msg: &'static str) {
            display("{} init error: {}", proj_name,  msg)
        }
        /// Encountered when mandatory keywords have not been found
        MandatoryWCSKeywordsMissing(keyword: &'static str) {
            display("{} keyword is mandatory for defining a WCS", keyword)
        }
        /// Not implemented projection error
        NotImplementedProjection(proj_name: String) {
            display("Projection '{}' has not been implemented!", proj_name)
        }
        /// Unrecognized radesys error
        UnrecognizedRadeSys(radesys: String) {
            display("Unrecognized radesys keyword value: '{}'!", radesys)
        }
        /// Celestial pole cannot be determined
        CelestialPoleInvalid {
            display("Bad combination of CRVAL1, CRVAL2, LONPOLE leading to an invalid celestial pole")
        }
    }
}
