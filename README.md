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
      * [Optional Pre-alignment](#optional-pre-alignment)
      * [Alignment Scoring](#alignment-scoring)
      * [Alignment mode](#alignment-mode)
      * [Circular contigs](#circular-contigs)
      * [Additional Output Options](#additional-output-options)
  
<!---toc end-->

## Disclaimer

**Stitch is under active development, and is currently alpha quality software for research purposes only.** Many features have not been tested, implemented, or considered. Please contact [Fulcrum Genomics](www.fulcrumgenomics.com) if you're considering using this software. 
Please submit an [issue](https://github.com/fulcrumgenomics/stitch/issues) - or better yet, a [pull request](https://github.com/fulcrumgenomics/stitch/pulls) - if you discover a bug or identify a missing feature.

## Overview

Stitch is a toolkit for analysis of chimeric reads in sequencing data (especially long reads like ONT and PacBio). 
Chimeric sequencing reads derive from template molecules that contain disjoint regions of the linear (or circular) reference genome, e.g., due to recombination, double-strand break repair, or gene editing. 
Chimeric reads may span large inter- or intra-contig distances, and thus present a challenge to traditional sequence aligners.

Potential use cases include, but are not limited to:

1. Aligning vectors/plasmids/viruses (assemblies, long reads, etc) to a database of constructs to determine the component structure
2. Determining the structure of complex structural variants (e.g., chromthripsis, cccDNA, etc.)
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
One record is marked as primary (see `--pick-primary`), while the remaining records are marked as secondary.  
The HI/HN SAM tags are used to denote the order in which the alignments occur.  
Furthermore, the records are output by this tool in the order in which they align to the query/read sequence.

Multiple contigs in the input FASTA are supported.

### Optional pre-alignment

The `-p`/`--pre-alignment` option may be used to pre-align the read using banded local alignment to select only those reads with a minimum pre-alignment score to perform the full alignment. 
Additional options are availabe to control the k-mer size, band-width, and minimum alignment score for this step.
Furthermore, the full alignmnet can be limited to align the read to just those reference contigs with minimum alignment score with the `-x`/`--pre-align-subset-contigs`.  
This is useful when aligning to a large database of individual constructs as contigs, where the read is expected only to align to a small subset of the contigs.

### Alignment Scoring

Alignments are scored using a match score, mismatch penalty, and affine gap penalty (a gap of size `k` costs `{-O} + {-E}*k`).  
Scores must be positive while penalities must be negative.

The jump score can be specified with `--jump-score`.

The jump score may also be specified specific to the jump being within
- the same contig and the same strand (`--jump-score-same-contig-and-strand`)
- the same contig but opposite strand (`--jump-score-same-contig-opposite-strand`), and
- the across different contigs (`--jump-score-inter-contig`).  

If any of these options are not specified, then they will default to the the value specified by `--jump-score`.

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

Additional options to pick which alignment is primary (versus secondary) for alignments that jump across contigs, as well as filtering out poorly-aligned secondary alignments are available 
