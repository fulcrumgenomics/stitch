use clap::Parser;
use env_logger::Env;
use fqcv_lib::util::version::built_info;

pub mod commands;

use anyhow::Result;
use commands::command::Command;
use enum_dispatch::enum_dispatch;
use fqcv_lib::commands::align::Align;

// #[global_allocator]
// static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[derive(Parser, Debug)]
struct Args {
    #[clap(subcommand)]
    subcommand: Subcommand,
}

#[enum_dispatch(Command)]
#[derive(Parser, Debug)]
#[command(version = built_info::VERSION.as_str())]
enum Subcommand {
    Align(Align),
}

fn main() -> Result<()> {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    let args: Args = Args::parse();
    args.subcommand.execute()
}
