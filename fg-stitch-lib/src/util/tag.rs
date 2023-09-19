//! Custom SAM tags for stitch.  See the README.md file for the description of each tag.
use noodles::sam::record::data::field::Tag;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum CustomTag {
    /// (`qs`)
    QueryStart,
    /// (`qe`)
    QueryEnd,
    /// (`ts`)
    TargetStart,
    /// (`te`)
    TargetEnd,
    /// (`as`)
    ChainAlignmentScore,
    /// (`xs`)
    SuboptimalScore,
    /// (`si`)
    SubAlignmentIndex,
    /// (`sc`)
    SubAlignmentCigar,
    /// (`cl`)
    ChainLength,
    /// (`ci`)
    ChainIndex,
    /// (`cn`)
    NumberOfChains,
}

impl AsRef<[u8]> for CustomTag {
    fn as_ref(&self) -> &[u8] {
        match self {
            CustomTag::QueryStart => b"qs",
            CustomTag::QueryEnd => b"qe",
            CustomTag::TargetStart => b"ts",
            CustomTag::TargetEnd => b"te",
            CustomTag::ChainAlignmentScore => b"as",
            CustomTag::SuboptimalScore => b"xs",
            CustomTag::SubAlignmentIndex => b"si",
            CustomTag::SubAlignmentCigar => b"sc",
            CustomTag::ChainLength => b"cl",
            CustomTag::ChainIndex => b"ci",
            CustomTag::NumberOfChains => b"cn",
        }
    }
}

impl From<CustomTag> for Tag {
    fn from(value: CustomTag) -> Self {
        value.as_ref().try_into().unwrap()
    }
}
