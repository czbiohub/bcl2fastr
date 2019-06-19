#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_xml_rs;

mod structs;
mod parser;


fn main() {
    let filename = "/mnt/bcl2fastr_test_data/190414_A00111_0296_AHJCWWDSXX/test_runinfo.xml".to_string();
    parser::main(filename);    
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
