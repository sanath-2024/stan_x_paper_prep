/// custom error type
#[derive(Debug)]
pub enum Error {
    /// CLI error
    CLI,
    /// file not found error
    FileNotFound,
    /// I/O error (besides file not found)
    IO,
    /// FASTA error
    FASTA,
    /// FASTQ error
    FASTQ,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for Error {}

/// custom result type
pub type Result<T> = std::result::Result<T, Error>;
