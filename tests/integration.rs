//! Integration testing for the CLI

#[cfg(test)]
mod integration {
    use assert_cli;

    #[test]
    fn run() {
        assert_cli::Assert::main_binary()
            .with_args(&["--run-path", "test_data/190414_A00111_0296_AHJCWWDSXX"])
            .with_args(&[
                "--samplesheet",
                "test_data/190414_A00111_0296_AHJCWWDSXX/SampleSheet.csv"
            ])
            .with_args(&["--output", "test_data/test_output"])
            .succeeds()
            .unwrap();
    }

    #[test]
    fn call_without_args() {
        assert_cli::Assert::main_binary()
            .fails()
            .and()
            .stderr()
            .contains("The following required arguments were not provided:")
            .unwrap();
    }
    
    #[test]
    fn no_run() {
        assert_cli::Assert::main_binary()
            .with_args(&["--run-path", "test_data/190414_A00111_0296_AHJCWWDSXXX"])
            .with_args(&[
                "--samplesheet",
                "test_data/190414_A00111_0296_AHJCWWDSXX/SampleSheet.csv"
            ])
            .with_args(&["--output", "test_data/test_output"])
            .fails()
            .and()
            .stderr()
            .contains("Could not find run path test_data/190414_A00111_0296_AHJCWWDSXXX")
            .unwrap();
    }

    #[test]
    fn no_samplesheet() {
        assert_cli::Assert::main_binary()
            .with_args(&["--run-path", "test_data/190414_A00111_0296_AHJCWWDSXX"])
            .with_args(&["--samplesheet", "test_data/no_file.csv"])
            .with_args(&["--output", "test_data/test_output"])
            .fails()
            .and()
            .stderr()
            .contains("Could not find samplesheet test_data/no_file.csv")
            .unwrap();
    }

    #[test]
    fn no_output() {
        assert_cli::Assert::main_binary()
            .with_args(&["--run-path", "test_data/190414_A00111_0296_AHJCWWDSXX"])
            .with_args(&[
                "--samplesheet",
                "test_data/190414_A00111_0296_AHJCWWDSXX/SampleSheet.csv"
            ])
            .with_args(&["--output", "test_data/not_a_dir"])
            .fails()
            .and()
            .stderr()
            .contains("Could not find output path test_data/not_a_dir")
            .unwrap();
    }

    #[test]
    fn bad_tile_chunk() {
        assert_cli::Assert::main_binary()
            .with_args(&["--run-path", "test_data/190414_A00111_0296_AHJCWWDSXX"])
            .with_args(&[
                "--samplesheet",
                "test_data/190414_A00111_0296_AHJCWWDSXX/SampleSheet.csv"
            ])
            .with_args(&["--output", "test_data/test_output"])
            .with_args(&["--tile-chunk", "not_a_number"])
            .fails()
            .and()
            .stderr()
            .contains("Invalid value: The argument 'not_a_number' isn't a valid value")
            .unwrap();
    }

}