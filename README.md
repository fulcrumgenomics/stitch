# stitch

<p align="center">
  <a href="https://github.com/fulcrumgenomics/stitch/actions?query=workflow%3ACheck"><img src="https://github.com/fulcrumgenomics/stitch/actions/workflows/build_and_test.yml/badge.svg" alt="Build Status"></a>
  <br>
</p>

<!---toc start-->
* [stitch](#stitch)
   * [Disclaimer](#disclaimer)
   * [Overview](#overview)
   * [Installing](#installing)
      * [Installing with cargo](#installing-with-cargo)
      * [Building From Source](#building-from-source)
   * [Usage](#usage)
      * [stitch align](#stitch-align)
      * [SAM Flags and Tags](#sam-flags-and-tags)
      * [Optional Pre-alignment](#optional-pre-alignment)
      * [Alignment Scoring](#alignment-scoring)
      * [Alignment mode](#alignment-mode)
      * [Circular contigs](#circular-contigs)
      * [Additional Output Options](#additional-output-options)
  
<!---toc end-->

<p>
<a href float="left"="https://fulcrumgenomics.com"><img src=".github/logos/fulcrumgenomics.svg" alt="Fulcrum Genomics" height="100"/></a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with stitch and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-brightgreen.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-blue.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

## Disclaimer

**Stitch is under active development, and is currently alpha quality software for research purposes only.** Many features have not been tested, implemented, or considered. Please contact [Fulcrum Genomics](www.fulcrumgenomics.com) if you're considering using this software.
Please submit an [issue](https://github.com/fulcrumgenomics/stitch/issues) - or better yet, a [pull request](https://github.com/fulcrumgenomics/stitch/pulls) - if you discover a bug or identify a missing feature.

## Overview

Stitch is a toolkit for analysis of chimeric reads in sequencing data (especially long reads like ONT and PacBio).
Chimeric sequencing reads derive from template molecules that contain disjoint regions of the linear (or circular) reference genome, e.g., due to recombination, double-strand break repair, or gene editing. 
Chimeric reads may span large inter- or intra-contig distances, and thus present a challenge to traditional sequence aligners.

Potential use cases include, but are not limited to:

1. Aligning vectors/plasmids/viruses (assemblies, long reads, etc) to a database of constructs to determine the component structure
2. Determining the structure of complex structural variants (e.g., chromothripsis, cccDNA, etc.)
3. Per-read evidence for somatic fusions and other structural variants

## Installing

### Installing with `cargo`
To install with cargo you must first [install rust](https://doc.rust-lang.org/cargo/getting-started/installation.html).

Then, to install `stitch` run:

```console
cargo install fg-stitch
```

### Building from source

First, clone the git repo:

```console
git clone https://github.com/fulcrumgenomics/stitch.git
```

Second, [install rust](https://doc.rust-lang.org/cargo/getting-started/installation.html) if you have not already.

Then build the toolkit in release mode:

```console
cd stitch
cargo build --release
./target/release/stitch --help
```

## Usage

### `stitch align`

Aligns reads of any length (in FASTQ format) against an assembly (in FASTA format).

Note that `stitch` was originally developed for the purpose of aligning reads from vectors/plasmids/viruses to a large number of relatively short constructs to determine each read's component structure.
The assembly is indexed at runtime and the index is held in memory, so performance and memory usage likely will not scale well to large genomes.

The alignment extends the traditional alignment algorithms by introducing a "jump" move/operator, whereby the alignment is able to jump anywhere in the reference sequence for a fixed cost.
This models the case where the read sequence is composed of multiple discrete segments of the expected reference, in shuffled, inverted, or truncated order.
In the default mode, the alignment is allowed to jump on the same strand either before or after the current reference position, whereas with `--double-strand` the alignment is also able to jump and continue on the opposite strand.
Since the alignment may jump to a previous reference position, different segments of the read may align to the same reference sequence multiple times.

The output is in SAM/BAM format.
If the alignment has `N` jumps, then the output will contain `N+1` records for the input read.
In some cases, multiple "chains" of alignments are output for a single read.
See [SAM Flags and Tags](#sam-flags-and-tags) for details on how the SAM flags are set as well as
custom SAM tags, to reason which alignments are part of which linear alignment chain, and which
linear alignment chain is the "primary".

Multiple contigs in the input FASTA are supported (and encouraged).

### SAM Flags and Tags

In the simplest case, a read (query) aligns to one contig (target) sub-sequence, such that a single
SAM alignment record is output.
If a read has a linear alignment (chain), whereby (typically) non-overlapping sub-sequences of the
read are aligned to various sub-sequences of the reference (may be different contigs), then one
SAM alignment record is output per sub-alignment in the chain.
In some cases, there may be multiple chains of alignments reported, such that one linear chain
is the "representative" (or best) linear alignment of the read, and the others are "secondary".
In these cases, SAM alignment records are output for each sub-alignment across each chain.

For example, there may be three (3) linear alignments of the read, with the first having six (6)
sub-alignments (i.e., five (5) jumps), the second having two (2) sub-alignments (i.e., one jump),
and the third having a single sub-alignment (i.e., no jumps).
In this case, there will be nine (9) SAM alignment records output (6 + 2 + 1).

SAM flags and custom SAM tags are used to identify the chain to which the SAM alignment record belongs, as 
well as the order of the alignment records in the given linear chain.  In particular, SAM flags
are set as follows:
1. One sub-alignment will (a) not have the secondary flag set, and (b) not have the
   supplementary flag set.  This is the "primary" (representative) sub-alignment in the
   "primary chain" (NB: not a formal term).
2. The remaining sub-alignments in the "primary chain" have the supplementary flag set,
   but not the secondary flag set.  They are part of the best linear alignment chain,
   but not the representative sub-alignment.
3. For "secondary chains", one sub-alignment is representative and has the secondary
   flag set (since it's not in the "primary chain") but not the supplementary flag.
4. For "secondary chains", the non-representative sub-alignments have both the
   secondary and supplementary flags set.

See `--pick-primary` for how the primary alignment is chosen in each chain.

The SAM tags are set as follows:

| tag | type | description |
| --- | --- | --- |
| `qs` | `i` | the zero-based index of the first query base in the sub-alignment |
| `qe` | `i` | the zero-based exclusive index of the last query base in the sub-alignment |
| `ts` | `i` |  the zero-based index of the first target base in the sub-alignment |
| `te` | `i` |  the zero-based exclusive index of the last target base in the sub-alignment |
| `as` | `i` |  the alignment score of the chain (not the sub-alignment, see `AS` for that) |
| `xs` | `i` |  the sub-optimal alignment score, practically the maximum of any pre-alignment and secondary chain |
| `si` | `i` |  the index of the sub-alignment in the current chain |
| `sc` | `Z` |  the cigar of the given sub-alignment, without any soft or hard clipping |
| `cl` | `i` |  the number of sub-alignments in the current chain |
| `ci` | `i` |  the index of the chain across all chains for this query |
| `cn` | `i` |  the number of chains for this query |
| `AS` | `i` |  the alignment score of the sub-alignment (not the chain, see `as` for that) |
| `SA` | `Z` |  the semicolon-delimited list of alignments for the given chain (rname, pos, strand, CIGAR, mapQ, NM) |
| `NM` | `i` | the number of edits in the sub-alignment |

The `-p`/`--pre-alignment` option may be used to pre-align the read using banded local alignment to select only those reads with a minimum pre-alignment score to perform the full alignment.
Additional options are available to control the k-mer size, band-width, and minimum alignment score for this step.
Furthermore, the full alignment can be limited to align the read to just those reference contigs with minimum alignment score with the `-x`/`--pre-align-subset-contigs`.
This is useful when aligning to a large database of individual constructs as contigs, where the read is expected only to align to a small subset of the contigs.

### Alignment Scoring

Alignments are scored using a match score, mismatch penalty, and affine gap penalty (a gap of size `k` costs `{-O} + {-E}*k`).
Scores must be positive while penalties must be negative.

The jump score can be specified with `--jump-score`.

The jump score may also be specified specific to the jump being within
- the same contig and the same strand (`--jump-score-same-contig-and-strand`)
- the same contig but opposite strand (`--jump-score-same-contig-opposite-strand`), and
- the across different contigs (`--jump-score-inter-contig`).

If any of these options are not specified, then they will default to the value specified by `--jump-score`.

### Alignment mode

Four major modes of alignment are supported with the `-m`/`--mode` option:

1. Local: aligns a sub-sequence of the read versus a sub-sequence of the reference.
2. QueryLocal: aligns a sub-sequence of the read versus the full reference.
3. TargetLocal: aligns the full read versus a sub-sequence of the reference.
4. Global: aligns the full read versus the full reference.

### Circular contigs

The `-C`/`--circular` option may be used to treat all input contigs as circular.
This allows the alignment to jump back to the start of the same contig at no cost (only from the end of the same contig).
In some alignment modes, a prefix or suffix of the read may be unaligned and is near the start or end of the contig.
In this case the read is re-aligned to better represent circular contigs, with `--circular-slop` controlling how close to the start/end of the contig the unaligned prefix/suffix must be for this to be triggered.

### Additional output options

Additional options to pick which alignment is primary (versus secondary) for alignments that jump across contigs, as well as filtering out poorly-aligned secondary alignments are available...
