//! Integration testing for the CLI

#[cfg(test)]
mod integration {
    use assert_cmd::crate_name;
    use assert_cmd::prelude::*;
    use predicates::prelude::*;
    use std::process::Command;

    #[test]
    fn run() {
        let mut cmd = Command::cargo_bin(crate_name!()).unwrap();
        cmd.args(&[
            "--run-path",
            "test_data/190414_A00111_0296_AHJCWWDSXX",
            "--samplesheet",
            "test_data/190414_A00111_0296_AHJCWWDSXX/SampleSheet.csv",
            "--output",
            "test_data/test_output",
        ]);

        cmd.assert().success();
    }

    #[test]
    fn call_without_args() {
        let mut cmd = Command::cargo_bin(crate_name!()).unwrap();

        cmd.assert().failure().stderr(
            predicate::str::contains("The following required arguments were not provided:")
                .from_utf8(),
        );
    }

    #[test]
    fn no_run() {
        let mut cmd = Command::cargo_bin(crate_name!()).unwrap();
        cmd.args(&[
            "--run-path",
            "test_data/190414_A00111_0296_AHJCWWDSXXX",
            "--samplesheet",
            "test_data/190414_A00111_0296_AHJCWWDSXX/SampleSheet.csv",
            "--output",
            "test_data/test_output",
        ]);

        cmd.assert().failure().stderr(
            predicate::str::contains(
                "Could not find run path test_data/190414_A00111_0296_AHJCWWDSXXX",
            )
            .from_utf8(),
        );
    }

    #[test]
    fn no_samplesheet() {
        let mut cmd = Command::cargo_bin(crate_name!()).unwrap();
        cmd.args(&[
            "--run-path",
            "test_data/190414_A00111_0296_AHJCWWDSXX",
            "--samplesheet",
            "test_data/no_file.csv",
            "--output",
            "test_data/test_output",
        ]);

        cmd.assert().failure().stderr(
            predicate::str::contains("Could not find samplesheet test_data/no_file.csv")
                .from_utf8(),
        );
    }

    #[test]
    fn no_output() {
        let mut cmd = Command::cargo_bin(crate_name!()).unwrap();
        cmd.args(&[
            "--run-path",
            "test_data/190414_A00111_0296_AHJCWWDSXX",
            "--samplesheet",
            "test_data/190414_A00111_0296_AHJCWWDSXX/SampleSheet.csv",
            "--output",
            "test_data/not_a_dir",
        ]);

        cmd.assert().failure().stderr(
            predicate::str::contains("Could not find output path test_data/not_a_dir").from_utf8(),
        );
    }

    #[test]
    fn bad_read_chunk() {
        let mut cmd = Command::cargo_bin(crate_name!()).unwrap();
        cmd.args(&[
            "--run-path",
            "test_data/190414_A00111_0296_AHJCWWDSXX",
            "--samplesheet",
            "test_data/190414_A00111_0296_AHJCWWDSXX/SampleSheet.csv",
            "--output",
            "test_data/test_output",
            "--read-chunks",
            "not_a_number",
        ]);

        cmd.assert().failure().stderr(
            predicate::str::contains(
                "Invalid value: The argument 'not_a_number' isn't a valid value",
            )
            .from_utf8(),
        );
    }
}
