# fqcv

<p align="center">
  <a href="https://github.com/fulcrumgenomics/fqcv/actions?query=workflow%3ACheck"><img src="https://github.com/fulcrumgenomics/fqcv/actions/workflows/build_and_test.yml/badge.svg" alt="Build Status"></a>
  <br>
</p>

Alignment for vector quality control

<!---toc start-->
* [fqcv align](#fqcv-align)
  * [Optional Pre-alignment](#optional-pre-alignment)
  * [Alignment Scoring](#alignment-scoring)
  * [Alignment mode](#alignment-mode)
  * [Circular contigs](#circular-contigs)
  * [Additional Output Options](#additional-output-options)
* [Installing](#installing)
  * [Installing with cargo](#installing-with-cargo)
  * [Building From Source](#building-from-source)
* [Developing](#developing)
* [Releasing a New Version](#releasing-a-new-version)
  * [Pre-requisites](#pre-requisites)
  * [Prior to Any Release](#prior-to-any-release)
  * [Semantic Versioning](#semantic-versioning)
  * [Major Release](#major-release)
  * [Minor and Patch Release](#minor-and-patch-release)
  * [Release Candidate](#release-candidate)
  
<!---toc end-->

## fqcv align

Perfoms alignment of a long reads against one or more reference/expected vector/plasmid/construct(s).

The alignment extends the traditional alignment algorithms by introducing a "jump" move/operator, 
whereby the alignment is able to jump anywhere in the reference sequence for a fixed cost.  This 
models the case where the read sequence is composed of multiple discrete segments of the expected 
reference, in shuffled, inverted, or truncated order.  In the default mode, the alignment is
allowed to jump on the same strand either before or after the current reference position, where as
with `--double-strand` the alignment is also able to jump and continue on the opposite
strand and continue on the opposite strand.  Since the alignment may jump to a previous reference
position, different segments of the read may align to the same reference sequence multiple times.

The output is in SAM/BAM format.  If the alignment has `N` jumps, then the output will contain
`N-1` records for the input read.  One record is marked as primary (see `--pick-primary`),
while the remaining records are marked as secondary.  The HI/HN SAM tags are used to denote
the order in which the alignments occur.  Furthermore, the order of the records output by this
tool are in the order in which they align the query/read sequence.

Multiple contigs in the input FASTA are supported.

### Optional Pre-alignment

The `-p`/`--pre-alignment` option may be used to pre-align the read using banded local alignment
to select only those reads with a minimum pre-alignment score to perform the full alignment.  Additional
options are availabe to control the k-mer size, band-width, and minimum alignment score for this step.
Furthermore, the full alignmnet can be limited to align the read to just those reference contigs
with minimum alignment score with the `-x`/`--pre-align-subset-contigs`.  This is useful when 
aligning to a large database of individual constructs as contigs, where the read is expected only
to align to a small subset of the contigs.

### Alignment Scoring

Alignments are scored using a match score, mismatch penaltiy, and affine gap penalty (a gap of size 
`k` costs `{-O} + {-E}*k`).  Scores must be positive while penalities must be negative.

The jump score can be specified with `--jump-score`.

The jump score may also be specified specific to the jump being within the same contig
and the same strand (`--jump-score-same-contig-and-strand`), the same contig but opposite
strand (`--jump-score-same-contig-opposite-strand`), and the across different contigs.
(`--jump-score-inter-contig`).  If any of these options are not specified, then they will
default to the the value specified by `--jump-score`.

### Alignment mode

Four major modes of alignment are supported with the `-m`/`--mode` option:

1. Local: aligns a sub-sequence of the read versus a sub-sequence of the reference.
2. QueryLocal: aligns a sub-sequence of the read versus the full reference.
3. TargetLocal: aligns the full read versus a sub-sequence of the reference.
4. Global: aligns the full read versus the full reference.

### Circular contigs

The `-C`/`--circular` option may be used to treat all input contigs as circular.  This allows the
alignment to jump back to the start of the same contig at no cost.  In some alignment modes, a
prefix or suffix of the read may be unaligned and is near the start or end of the contig.  In this
case the read is re-aligned to better represent circular contigs, with `--circular-slop` controlling
how near the unaligned prefix/suffix must be to the start/end of the contig for this to be triggered.

### Additional Output Options

Additional options to pick which alignment is primary (versus secondary) for alignments that jump
across contigs, as well as filtering out poorly aligned secondary alignments are given.

## Installing

### Installing with `cargo`
To install with cargo you must first [install rust](https://doc.rust-lang.org/cargo/getting-started/installation.html).
Which (On Mac OS and Linux) can be done with the command:

```console
curl https://sh.rustup.rs -sSf | sh
```

Then, to install `fqcv` run:

```console
cargo install fqcv
```

### Building From Source

First, clone the git repo:

```console
git clone https://github.com/fulcrumgenomics/fqcv.git
```

Secondly, if you do not already have rust development tools installed, install via [rustup](https://rustup.rs/):

```console
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Then build the toolkit in release mode:

```console
cd fqcv
cargo build --release
./target/release/fqcv --help
```

## Developing

fqcv is developed in Rust and follows the conventions of using `rustfmt` and `clippy` to ensure both code quality and standardized formatting.
When working on fqcv, before pushing any commits, please first run `./ci/check.sh` and resolve any issues that are reported.

## Releasing a New Version

### Pre-requisites

Install [`cargo-release`][cargo-release-link]

```console
cargo install cargo-release
```

### Prior to Any Release

Create a release that will not try to push to `crates.io` and verify the command:

```console
cargo release [major,minor,patch,release,rc...] --no-publish
```

Note: "dry-run" is the default for cargo release.

See the [`cargo-release` reference documentation][cargo-release-docs-link] for more information

### Semantic Versioning

This tool follows [Semantic Versioning](https://semver.org/).  In brief:

* MAJOR version when you make incompatible API changes,
* MINOR version when you add functionality in a backwards compatible manner, and
* PATCH version when you make backwards compatible bug fixes.

### Major Release

To create a major release:

```console
cargo release major --execute
```

This will remove any pre-release extension, create a new tag and push it to github, and push the release to creates.io.

Upon success, move the version to the [next candidate release](#release-candidate).

Finally, make sure to [create a new release][new-release-link] on GitHub.

### Minor and Patch Release

To create a _minor_ (_patch_) release, follow the [Major Release](#major-release) instructions substituting `major` with `minor` (`patch`):

```console
cargo release minor --execute
```

### Release Candidate

To move to the next release candidate:

```console
cargo release rc --no-tag --no-publish --execute
```

This will create or bump the pre-release version and push the changes to the main branch on github.
This will not tag and publish the release candidate.
If you would like to tag the release candidate on github, remove `--no-tag` to create a new tag and push it to github.

[cargo-release-link]:      https://github.com/crate-ci/cargo-release
[cargo-release-docs-link]: https://github.com/crate-ci/cargo-release/blob/master/docs/reference.md
[new-release-link]:        https://github.com/fulcrumgenomics/fqcv/releases/new
