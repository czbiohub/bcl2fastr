#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_xml_rs;

mod parser;


fn main() {
    let filename_info = "src/test_data/RunInfo.xml".to_string();
    let filename_params = "src/test_data/test_runparams.xml".to_string();
    parser::parse_run_info(filename_info);
    parser::parse_run_params(filename_params);
}


fn version_info() -> &'static str {
    "bcl2fastr beta version"
}




#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_version_info() {
        assert_eq!(version_info(), "bcl2fastr beta version");
    }
}
