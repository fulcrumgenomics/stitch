use anyhow::Result;
use clap::builder::PossibleValue;
use enum_dispatch::enum_dispatch;
use std::{fmt::Display, str::FromStr};

#[enum_dispatch]
pub trait Command {
    #[allow(clippy::missing_errors_doc)]
    fn execute(&self) -> Result<()>;
}

pub trait ValueEnum: Display + FromStr {
    fn variants<'a>() -> &'a [Self];

    fn possible_values() -> Vec<PossibleValue> {
        Self::variants()
            .iter()
            .map(|variant| PossibleValue::new(variant.to_string()))
            .collect()
    }
}
