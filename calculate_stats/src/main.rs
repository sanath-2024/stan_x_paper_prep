//! calculate_stats: calculate various statistics about our reads and the genome in general

#![warn(missing_docs)]

mod cli;
mod fast_coverage;
mod result;
mod util;

pub use result::{Error, Result};

use cli::Command;

fn main() -> Result<()> {
    let matches = cli::get().get_matches();
    let (name, matches) = matches.subcommand().unwrap();
    match name {
        "fast_coverage" => {
            fast_coverage::FastCoverage::run(fast_coverage::FastCoverage::read_matches(matches)?)
        }
        _ => Err(Error::CLI),
    }
}
