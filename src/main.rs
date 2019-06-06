fn main() {

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
