use clap::{Arg, ArgMatches, ValueHint};

/// trait Command represents a set of parsed arguments that can be used to run a command
pub trait Command: Sized {
    fn read_matches(matches: &ArgMatches) -> crate::Result<Self>;
    fn run(args: Self) -> crate::Result<()>;
}

/// reusable argument for the fastq reads
fn reads_arg() -> Arg {
    Arg::new("reads")
        .long("reads")
        .short('q')
        .help("path to FASTQ reads file (\"-\" means standard input)")
        .value_hint(ValueHint::FilePath)
}

/// reusable argument for the reference genome
fn ref_arg() -> Arg {
    Arg::new("ref")
        .long("ref")
        .short('r')
        .help("path to FASTA reference file (\"-\" means standard input)")
        .value_hint(ValueHint::FilePath)
}

/// reusable argument for an output file
fn output_arg() -> Arg {
    Arg::new("out")
        .long("out")
        .short('o')
        .help("path to output file (\"-\" means standard output)")
        .value_hint(ValueHint::FilePath)
}

pub fn get() -> clap::Command {
    clap::Command::new("Stan-X: Calculate Statistics")
        .author("Sanath Govindarajan <sgovindarajan@utexas.edu>")
        .version(clap::crate_version!())
        .about("calculate various statistics for the Stan-X paper in a reproducible matter")
        .subcommand(
            clap::Command::new("fast_coverage")
                .about("calculate coverage in the naive way (number of nt in reads / number of nt in genome)")
                .arg(reads_arg().required(true))
                .arg(ref_arg().required(true))
                .arg(output_arg().required(true))
        )
        .subcommand_required(true)
}
