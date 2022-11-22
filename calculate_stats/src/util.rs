use crate::{Error, Result};

use std::{
    fs::File,
    io::{stdin, stdout, Read, Stdin, Stdout, Write},
    path::Path,
};

pub enum InFile {
    Regular(File),
    Stdin(Stdin),
}

impl Read for InFile {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            InFile::Regular(f) => f.read(buf),
            InFile::Stdin(s) => s.read(buf),
        }
    }
}

pub enum OutFile {
    Regular(File),
    Stdout(Stdout),
}

impl Write for OutFile {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            OutFile::Regular(f) => f.write(buf),
            OutFile::Stdout(s) => s.write(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            OutFile::Regular(f) => f.flush(),
            OutFile::Stdout(s) => s.flush(),
        }
    }
}

pub fn open_file_read<P: AsRef<Path>>(path: P) -> Result<InFile> {
    if path.as_ref().to_str().unwrap() == "-" {
        return Ok(InFile::Stdin(stdin()));
    }
    Ok(InFile::Regular(File::open(path).map_err(
        |e| match e.kind() {
            std::io::ErrorKind::NotFound => Error::FileNotFound,
            _ => Error::IO,
        },
    )?))
}

pub fn create_file_write<P: AsRef<Path>>(path: P) -> Result<OutFile> {
    if path.as_ref().to_str().unwrap() == "-" {
        return Ok(OutFile::Stdout(stdout()));
    }
    Ok(OutFile::Regular(File::create(path).map_err(
        |e| match e.kind() {
            std::io::ErrorKind::NotFound => Error::FileNotFound,
            _ => Error::IO,
        },
    )?))
}
