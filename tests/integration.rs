//! Integration testing for the CLI

#[cfg(test)]
mod integration {
    use assert_cli;

    #[test]
    fn call_without_args() {
        assert_cli::Assert::main_binary()
            .fails()
            .and()
            .stderr()
            .contains("error: The following required arguments were not provided:")
            .unwrap();
    }
    
    #[test]
    fn no_run_() {
        assert_cli::Assert::main_binary()
            .with_args(&["--run-path", "test_data/190414_A00111_0296_AHJCWWDSXXX"])
            .with_args(&["--samplesheet", "test_data/190414_A00111_0296_AHJCWWDSXX/190414_A00111_0296_AHJCWWDSXX.csv"])
            .with_args(&["--output", "test_data/test_output"])
            .fails()
            .and()
            .stderr()
            .contains("Error reading NovaSeq run: No such file or directory")
            .unwrap();
    }
}