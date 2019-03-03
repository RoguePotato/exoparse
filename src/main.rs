//==============================================================================
//! # exoparse
//!
//! `exoparse` is a tool to parse exoplanet data.
//==============================================================================

extern crate exoparse;

use std::env;
use std::process;

use exoparse::Config;

fn main() {
    // Create a configuration based on the command line parameters provided.
    let config = Config::new(env::args()).unwrap_or_else(|err| {
        eprintln!("Problem parsing arguments: {}", err);
        process::exit(1);
    });

    // Run the program with the current configuration. If an error is returned
    // from run, then report it and exit.
    if let Err(err) = exoparse::run(config) {
        eprintln!("Application error: {}", err);
        process::exit(1);
    }
}
