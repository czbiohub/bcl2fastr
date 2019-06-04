fn main() {
    println!("Hello, world!");
}



#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_example() {
        assert_eq!(main(), "Hello ");
    }
}
