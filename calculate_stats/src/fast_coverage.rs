use crate::{util, Error, Result};

use rayon::prelude::*;

use std::io::Write;

pub struct FastCoverage {
    reads: util::InFile,
    reference: util::InFile,
    out: util::OutFile,
}

impl crate::cli::Command for FastCoverage {
    fn read_matches(matches: &clap::ArgMatches) -> Result<Self> {
        let reads = matches.get_one::<String>("reads").ok_or(Error::CLI)?;
        let reference = matches.get_one::<String>("ref").ok_or(Error::CLI)?;
        let out = matches.get_one::<String>("out").ok_or(Error::CLI)?;
        let reads = util::open_file_read(reads)?;
        let reference = util::open_file_read(reference)?;
        let out = util::create_file_write(out)?;
        Ok(FastCoverage {
            reads,
            reference,
            out,
        })
    }
    fn run(mut args: Self) -> Result<()> {
        let reference = bio::io::fasta::Reader::new(args.reference);
        let ref_nt = reference
            .records()
            .par_bridge()
            .map(|rec| {
                let rec = rec.map_err(|_| Error::FASTQ)?;
                Ok(rec.seq().len())
            })
            .reduce(
                || Ok(0),
                |a: Result<usize>, b| {
                    let a = a?;
                    let b = b?;
                    Ok(a + b)
                },
            )?;
        write!(args.out, "reference nt: {}\n", ref_nt).map_err(|_| Error::IO)?;
        let reads = bio::io::fastq::Reader::new(args.reads);
        let reads_nt = reads
            .records()
            .par_bridge()
            .map(|rec| {
                let rec = rec.map_err(|_| Error::FASTQ)?;
                Ok(rec.seq().len())
            })
            .reduce(
                || Ok(0),
                |a: Result<usize>, b| {
                    let a = a?;
                    let b = b?;
                    Ok(a + b)
                },
            )?;
        write!(args.out, "reads nt: {}\n", reads_nt).map_err(|_| Error::IO)?;
        write!(args.out, "coverage: {}x\n", reads_nt as f64 / ref_nt as f64)
            .map_err(|_| Error::IO)?;
        Ok(())
    }
}
