use fqcv_lib::opts::setup;
use fqcv_lib::run::run;
use log::error;
use std::process::exit;

fn main() {
    let opts = setup();

    if let Err(err) = run(&opts) {
        error!("{:#}", err);
        exit(1);
    }
}
