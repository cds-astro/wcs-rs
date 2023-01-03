use fitsrs::hdu::Header;

use crate::error::Error;

pub unsafe fn string_to_keyword_type(keyword: &str) -> &[u8; 8] {
    keyword.as_bytes()[0..=7]
        .try_into()
        .expect("String has less than 8 characters")
}

pub fn retrieve_mandatory_parsed_keyword<T>(header: &Header, keyword: &'static str) -> Result<T, Error>
where
    T: fitsrs::card::CardValue
{
    if let Some(val) = header.get_parsed::<T>(unsafe { string_to_keyword_type(keyword) } ) {
        val.map_err(|e| e.into())
    } else {
        Err(Error::MandatoryWCSKeywordsMissing(keyword))
    }
}