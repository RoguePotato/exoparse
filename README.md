# Exoparse

Parses the CSV formatted [EU](exoplanet.eu) or
[NASA](https://exoplanetarchive.ipac.caltech.edu/) exoplanet archives. Output is
provided as a consisent number of column, tab-separated values where NULL values
or string in the original database are replaced with zeros. This makes it easier
for plotting packages to read.

## Requirements

The Rust compiler, obtained via

```
$ curl https://sh.rustup.rs -sSf | sh
```

## Building and running

Building is performed with
```
cargo build
```
and running with
```
cargo run <database_name>
```