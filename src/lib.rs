//==============================================================================
//! # exoparse
//!
//! A containment of library functions which are required by exoparse
//==============================================================================

use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::prelude::*;

/// Definition for the configuration. The filename is the raw database file to
/// be read, the format (tab or CSV) defines how the data are presented. The
/// database is from which database is being read.
pub struct Config {
    pub filename: String,
    // TODO: Database, delimiter, included errors.
}

impl Config {
    /// Create a new configuration based on command line parameters. A result is
    /// returned based on whether the correct number of parameters are provided.
    pub fn new(mut args: std::env::Args) -> Result<Config, &'static str> {
        // Skip the first argument, this is the executable name.
        args.next();

        // Set the configuration parameters if they exist.
        let filename = match args.next() {
            Some(arg) => arg,
            None => return Err("No filename provided!"),
        };

        // Return the configuration with the chosen parameters.
        Ok(Config { filename })
    }
}

// Runs the main program, and (if they show up) catches errors along the way.
pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    // Read in the database.
    let contents = fs::read_to_string(config.filename)?;

    let mut data: Vec<Record> = Vec::new();

    // Go through all the lines of the data file.
    for line in contents.lines() {
        // Don't include commented lines or the header to columns.
        if line.chars().next().unwrap() != '#' && line.chars().next().unwrap() != 'r' {
            // Split the line by commas and somehow add it to the records.
            let mut split: Vec<&str> = line.split(",").collect();
            data.push(parse_line(split));
        }
    }

    // Fill in the NULL quantities.
    for i in 0..data.len() {
        for j in 0..data[i].fields.len() {
            if data[i].fields[j] == "" {
                data[i].fields[j] = "0".to_string();
            }
        }
    }

    // Write out the new file.
    let mut file = File::create("./output.dat")?;
    for i in 0..data.len() {
        for j in 0..data[i].fields.len() {
            if data[i].fields[j] != "^^^" {
                match data[i].fields[j].parse::<f32>() {
                    Ok(_ok) => {
                        file.write_fmt(format_args!("{}\t", data[i].fields[j]))?;
                    }
                    Err(_e) => {
                        file.write_all(b"0 ")?;
                    }
                }
            }
        }
        file.write_all(b"\n")?;
    }

    Ok(())
}

pub struct Record {
    fields: Vec<String>,
}

pub fn parse_line(line: Vec<&str>) -> Record {
    let mut fields: Vec<String> = vec![String::new(); line.len()];
    let mut skip: bool = false;

    // Go over each field and add it to the record. If the field begins with a
    // quote mark, start skipping until the closing quote is found. These quote
    // entries include extra commas. For the skipped entries, put a string ^^^
    // to be filtered out later.
    for i in 0..line.len() {
        if line[i].chars().count() > 0 {
            if line[i].chars().next().unwrap() == '"' {
                skip = true;
            }
            if line[i].chars().last().unwrap() == '"' {
                skip = false;
            }
            if !skip {
                fields[i] = line[i].to_string();
            } else {
                fields[i] = "^^^".to_string();
            }
        }
    }

    Record { fields }
}

// Exoplanets EU
// 000 - name
// 001 - planet_status
// 002 - mass
// 003 - mass_error_min
// 004 - mass_error_max
// 005 - mass_sini
// 006 - mass_sini_error_min
// 007 - mass_sini_error_max
// 008 - radius
// 009 - radius_error_min
// 010 - radius_error_max
// 011 - orbital_period
// 012 - orbital_period_error_min
// 013 - orbital_period_error_max
// 014 - semi_major_axis
// 015 - semi_major_axis_error_min
// 016 - semi_major_axis_error_max
// 017 - eccentricity
// 018 - eccentricity_error_min
// 019 - eccentricity_error_max
// 020 - inclination
// 021 - inclination_error_min
// 022 - inclination_error_max
// 023 - angular_distance
// 024 - discovered
// 025 - updated
// 026 - omega
// 027 - omega_error_min
// 028 - omega_error_max
// 029 - tperi
// 030 - tperi_error_min
// 031 - tperi_error_max
// 032 - tconj
// 033 - tconj_error_min
// 034 - tconj_error_max
// 035 - tzero_tr
// 036 - tzero_tr_error_min
// 037 - tzero_tr_error_max
// 038 - tzero_tr_sec
// 039 - tzero_tr_sec_error_min
// 040 - tzero_tr_sec_error_max
// 041 - lambda_angle
// 042 - lambda_angle_error_min
// 043 - lambda_angle_error_max
// 044 - impact_parameter
// 045 - impact_parameter_error_min
// 046 - impact_parameter_error_max
// 047 - tzero_vr
// 048 - tzero_vr_error_min
// 049 - tzero_vr_error_max
// 050 - k
// 051 - k_error_min
// 052 - k_error_max
// 053 - temp_calculated
// 054 - temp_calculated_error_min
// 055 - temp_calculated_error_max
// 056 - temp_measured
// 057 - hot_point_lon
// 058 - geometric_albedo
// 059 - geometric_albedo_error_min
// 060 - geometric_albedo_error_max
// 061 - log_g
// 062 - publication
// 063 - detection_type
// 064 - mass_detection_type
// 065 - radius_detection_type
// 066 - alternate_names
// 067 - molecules
// 068 - star_name
// 069 - ra
// 070 - dec
// 071 - mag_v
// 072 - mag_i
// 073 - mag_j
// 074 - mag_h
// 075 - mag_k
// 076 - star_distance
// 077 - star_distance_error_min
// 078 - star_distance_error_max
// 079 - star_metallicity
// 080 - star_metallicity_error_min
// 081 - star_metallicity_error_max
// 082 - star_mass
// 083 - star_mass_error_min
// 084 - star_mass_error_max
// 085 - star_radius
// 086 - star_radius_error_min
// 087 - star_radius_error_max
// 088 - star_sp_type
// 089 - star_age
// 090 - star_age_error_min
// 091 - star_age_error_max
// 092 - star_teff
// 093 - star_teff_error_min
// 094 - star_teff_error_max
// 095 - star_detected_disc
// 096 - star_magnetic_field
// 097 - star_alternate_names

// Structure of the Caltech exoplanet database with all fields.
// pub struct CaltechRecord {
// id: i32,                         // 000 - Database ID
// pl_hostname: Option<String>,     // 001 - Host Name
// pl_letter: Option<char>,         // 002 - Planet Letter
// pl_name: Option<String>,         // 003 - Planet Name
// pl_discmethod: Option<String>,   // 004 - Discovery Method
// pl_pnum: Option<u32>,            // 005 - Number of Planets in System
// pl_orbper: Option<f32>,          // 006 - Orbital Period [days]
// pl_orbsmax: Option<f32>,         // 007 - Orbit Semi-Major Axis [AU])
// pl_orbeccen: Option<f32>,        // 008 - Eccentricity
// pl_orbincl: Option<f32>,         // 009 - Inclination [deg]
// pl_bmassj: Option<f32>,          // 010 - Planet Mass or M*sin(i) [Jupiter mass]
// pl_bmassprov: Option<String>,    // 011 - Planet Mass or M*sin(i) Provenance
// pl_radj: Option<f32>,            // 012 - Planet Radius [Jupiter radii]
// pl_dens: Option<f32>,            // 013 - Planet Density [g/cm**3]
// pl_ttvflag: Option<bool>,        // 014 - TTV Flag
// pl_kepflag: Option<bool>,        // 015 - Kepler Field Flag
// pl_k2flag: Option<bool>,         // 016 - K2 Mission Flag
// pl_nnotes: Option<u32>,          // 017 - Number of Notes
// ra_str: Option<String>,          // 018 - RA [sexagesimal]
// ra: Option<f32>,                 // 019 - RA [decimal degrees]
// dec_str: Option<String>,         // 020 - Dec [sexagesimal]
// dec: Option<f32>,                // 021 - Dec [decimal degrees]
// st_dist: Option<f32>,            // 022 - Distance [pc]
// st_optmag: Option<f32>,          // 023 - Optical Magnitude [mag]
// st_optband: Option<String>,      // 024 - Optical Magnitude Band
// gaia_gmag: Option<f32>,          // 025 - G-band (Gaia) [mag]
// st_teff: Option<f32>,            // 026 - Effective Temperature [K]
// st_mass: Option<f32>,            // 027 - Stellar Mass [Solar mass]
// st_rad: Option<f32>,             // 028 - Stellar Radius [Solar radii]
// rowupdate: Option<String>,       // 029 - Date of Last Update
// pl_tranflag: Option<bool>,       // 030 - Planet Transit Flag
// pl_imgflag: Option<bool>,        // 031 - Planet Imaging Flag
// pl_rvflag: Option<bool>,         // 032 - Planet RV Flag
// pl_astflag: Option<bool>,        // 033 - Planet Astrometry Flag
// pl_omflag: Option<bool>,         // 034 - Planet Orbital Modulation Flag
// pl_cbflag: Option<bool>,         // 035 - Planet Circumbinary Flag
// pl_angsep: Option<f32>,          // 036 - Calculated Angular Separation [mas]
// pl_orbtper: Option<f32>,         // 037 - Time of Periastron [days]
// pl_orblper: Option<f32>,         // 038 - Long. of Periastron [deg]
// pl_rvamp: Option<f32>,           // 039 - Radial Velocity Amplitude [m/s]
// pl_eqt: Option<f32>,             // 040 - Equilibrium Temperature [K]
// pl_insol: Option<f32>,           // 041 - Insolation Flux [Earth flux]
// pl_massj: Option<f32>,           // 042 - Planet Mass [Jupiter mass]
// pl_msinij: Option<f32>,          // 043 - Planet M*sin(i) [Jupiter mass]
// pl_masse: Option<f32>,           // 044 - Planet Mass [Earth mass]
// pl_msinie: Option<f32>,          // 045 - Planet M*sin(i) [Earth mass]
// pl_bmasse: Option<f32>,          // 046 - Planet Mass or M*sin(i) [Earth mass]
// pl_rade: Option<f32>,            // 047 - Planet Radius [Earth radii]
// pl_rads: Option<f32>,            // 048 - Planet Radius [Solar radii]
// pl_trandep: Option<f32>,         // 049 - Transit Depth [percent]
// pl_trandur: Option<f32>,         // 050 - Transit Duration [days]
// pl_tranmid: Option<f32>,         // 051 - Transit Midpoint [days]
// pl_tsystemref: Option<String>,   // 052 - Time System Reference
// pl_imppar: Option<f32>,          // 053 - Impact Parameter
// pl_occdep: Option<f32>,          // 054 - Occultation Depth [percentage]
// pl_ratdor: Option<f32>,          // 055 - Ratio of Distance to Stellar Radius
// pl_ratror: Option<f32>,          // 056 - Ratio of Planet to Stellar Radius
// pl_def_reflink: Option<String>,  // 057 - Default Reference
// pl_disc: Option<u32>,            // 058 - Year of Discovery
// pl_disc_reflink: Option<String>, // 059 - Discovery Reference
// pl_locale: Option<String>,       // 060 - Discovery Locale
// pl_facility: Option<String>,     // 061 - Discovery Facility
// pl_telescope: Option<String>,    // 062 - Discovery Telescope
// pl_instrument: Option<String>,   // 063 - Discovery Instrument
// pl_status: Option<u32>,          // 064 - Status
// pl_mnum: Option<u32>,            // 065 - Number of Moons in System
// pl_st_npar: Option<u32>,         // 066 - Number of Stellar and Planet Parameters
// pl_st_nref: Option<u32>,         // 067 - Number of Stellar and Planet References
// pl_pelink: Option<String>,       // 068 - Link to Exoplanet Encyclopaedia
// pl_edelink: Option<String>,      // 069 - Link to Exoplanet Data Explorer
// pl_publ_date: Option<String>,    // 070 - Publication Date
// hd_name: Option<String>,         // 071 - HD Name
// hip_name: Option<String>,        // 072 - HIP Name
// st_rah: Option<f32>,             // 073 - RA [hrs]
// st_glon: Option<f32>,            // 074 - Galactic Longitude [deg]
// st_glat: Option<f32>,            // 075 - Galactic Latitude [deg]
// st_elon: Option<f32>,            // 076 - Ecliptic Longitude [deg]
// st_elat: Option<f32>,            // 077 - Ecliptic Latitude [deg]
// st_plx: Option<f32>,             // 078 - Parallax [mas]
// gaia_plx: Option<f32>,           // 079 - Gaia Parallax [mas]
// gaia_dist: Option<f32>,          // 080 - Gaia Distance [pc]
// st_pmra: Option<f32>,            // 081 - Proper Motion (RA) [mas/yr]
// st_pmdec: Option<f32>,           // 082 - Proper Motion (Dec) [mas/yr]
// st_pm: Option<f32>,              // 083 - Total Proper Motion [mas/yr]
// gaia_pmra: Option<f32>,          // 084 - Gaia Proper Motion (RA) [mas/yr]
// gaia_pmdec: Option<f32>,         // 085 - Gaia Proper Motion (Dec) [mas/yr]
// gaia_pm: Option<f32>,            // 086 - Gaia Total Proper Motion [mas/yr]
// st_radv: Option<f32>,            // 087 - Radial Velocity [km/s]
// st_sp: Option<String>,           // 088 - Spectral Type TODO: Check this one.
// st_spstr: Option<String>,        // 089 - Spectral Type
// st_logg: Option<f32>,            // 090 - Stellar Surface Gravity [log10(cm/s**2)]
// st_lum: Option<f32>,             // 091 - Stellar Luminosity [log(Solar)]
// st_dens: Option<f32>,            // 092 - Stellar Density [g/cm**3]
// st_metfe: Option<f32>,           // 093 - Stellar Metallicity [dex]
// st_metratio: Option<String>,     // 094 - Metallicity Ratio [Fe/H] etc.
// st_age: Option<f32>,             // 095 - Stellar Age [Gyr]
// st_vsini: Option<f32>,           // 096 - Rot. Velocity V*sin(i) [km/s]
// st_acts: Option<f32>,            // 097 - Stellar Activity S-index
// st_actr: Option<f32>,            // 098 - Stellar Activity log(R'HK)
// st_actlx: Option<f32>,           // 099 - X-ray Activity log(L<sub>x</sub>)
// swasp_id: Option<String>,        // 100 - SWASP Identifier
// st_nts: Option<u32>,             // 101 - Number of Time Series
// st_nplc: Option<u32>,            // 102 - Number of Planet Transit Light Curves
// st_nglc: Option<u32>,            // 103 - Number of General Light Curves
// st_nrvc: Option<u32>,            // 104 - Number of Radial Velocity Time Series
// st_naxa: Option<u32>,            // 105 - Number of Amateur Light Curves
// st_nimg: Option<u32>,            // 106 - Number of Images
// st_nspec: Option<u32>,           // 107 - Number of Spectra
// st_uj: Option<f32>,              // 108 - U-band (Johnson) [mag]
// st_vj: Option<f32>,              // 109 - V-band (Johnson) [mag]
// st_bj: Option<f32>,              // 110 - B-band (Johnson) [mag]
// st_rc: Option<f32>,              // 111 - R-band (Cousins) [mag]
// st_ic: Option<f32>,              // 112 - I-band (Cousins) [mag]
// st_j: Option<f32>,               // 113 - J-band (2MASS) [mag]
// st_h: Option<f32>,               // 114 - H-band (2MASS) [mag]
// st_k: Option<f32>,               // 115 - Ks-band (2MASS) [mag]
// st_wise1: Option<f32>,           // 116 - WISE 3.4um [mag]
// st_wise2: Option<f32>,           // 117 - WISE 4.6um [mag]
// st_wise3: Option<f32>,           // 118 - WISE 12.um [mag]
// st_wise4: Option<f32>,           // 119 - WISE 22.um [mag]
// st_irac1: Option<f32>,           // 120 - IRAC 3.6um [mag]
// st_irac2: Option<f32>,           // 121 - IRAC 4.5um [mag]
// st_irac3: Option<f32>,           // 122 - IRAC 5.8um [mag]
// st_irac4: Option<f32>,           // 123 - IRAC 8.0um [mag]
// st_mips1: Option<f32>,           // 124 - MIPS 24um [mag]
// st_mips2: Option<f32>,           // 125 - MIPS 70um [mag]
// st_mips3: Option<f32>,           // 126 - MIPS 160um [mag]
// st_iras1: Option<f32>,           // 127 - IRAS 12um Flux [Jy]
// st_iras2: Option<f32>,           // 128 - IRAS 25um Flux [Jy]
// st_iras3: Option<f32>,           // 129 - IRAS 60um Flux [Jy]
// st_iras4: Option<f32>,           // 130 - IRAS 100um Flux [Jy]
// st_photn: Option<u32>,           // 131 - Number of Photometry Measurements
// st_umbj: Option<f32>,            // 132 - U-B (Johnson) [mag]
// st_bmvj: Option<f32>,            // 133 - B-V (Johnson) [mag]
// st_vjmic: Option<f32>,           // 134 - V-I (Johnson-Cousins) [mag]
// st_vjmrc: Option<f32>,           // 135 - V-R (Johnson-Cousins) [mag]
// st_jmh2: Option<f32>,            // 136 - J-H (2MASS) [mag]
// st_hmk2: Option<f32>,            // 137 - H-Ks (2MASS) [mag]
// st_jmk2: Option<f32>,            // 138 - J-Ks (2MASS) [mag]
// st_bmy: Option<f32>,             // 139 - b-y (Stromgren) [mag]
// st_m1: Option<f32>,              // 140 - m1 (Stromgren) [mag]
// st_c1: Option<f32>,              // 141 - c1 (Stromgren) [mag]
// st_colorn: Option<u32>,          // 142 - Number of Color Measurements
// }