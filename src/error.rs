quick_error! {
    #[derive(Debug)]
    pub enum Error {
        /// Encountered when parsing the Fits file, returned from fitsrs
        FitsParsing(err: fitsrs::error::Error) {
            from()
        }
        /// Encountered when the projection cannot be instantiated
        /// i.e. Zpn needs coefficients meeting a specific condition
        InitProjection(proj_name: &'static str, msg: &'static str) {
            display("{} init error: {}", proj_name,  msg)
        }
        /// Not implemented projection error
        NotImplementedProjection(proj_name: String) {
            display("Projection '{}' has not been implemented!", proj_name)
        }
    }
}
