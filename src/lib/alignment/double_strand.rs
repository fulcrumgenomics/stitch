// Original Code was copied from:
//     https://github.com/rust-bio/rust-bio/blob/master/src/alignment/pairwise/mod.rs
// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp::max;
use std::i32;
use std::iter::repeat;

use crate::alignment::constants::AlignmentMode;
use bio::utils::TextSlice;
use serde::Deserialize;
use serde::Serialize;

use crate::alignment::constants::AlignmentOperation;
use crate::alignment::constants::DEFAULT_ALIGNER_CAPACITY;
use crate::alignment::pairwise::PairwiseAlignment;
use crate::alignment::scoring::MatchFunc;
use crate::alignment::scoring::Scoring;
use crate::alignment::single_strand::SingleStrandMatrix;
use crate::alignment::single_strand::XBuffer;

use super::constants::MIN_SCORE;
use super::traceback::TB_DEL;
use super::traceback::TB_INS;
use super::traceback::TB_MATCH;
use super::traceback::TB_START;
use super::traceback::TB_SUBST;
use super::traceback::TB_XCLIP_PREFIX;
use super::traceback::TB_XCLIP_SUFFIX;
use super::traceback::TB_YCLIP_PREFIX;
use super::traceback::TB_YCLIP_SUFFIX;

/// A generalized Smith-Waterman aligner, allowing for the alignment to jump forward
/// (or backward) in `x` (on either same strand).
///
/// `M(i,j,k)` is the best score such that the alignment jumps from where `x[k]` and `y[j]` ends
/// in a match (or substitution), to where `x[i]` and `y[j]` ends in a match (or substitution).
/// `M(i,j,i-1)` is the equivalent to Smith-Waterman, whereas `M(i,j,k)` where `k != i - 1`
/// includes a jump from a position in `x` (can be before or after `i - 1`).
/// ```ignore
///              .... A   G  x_i
///              .... C   G  y_j
/// ```
///
/// `I(i,j)` is the best score such that `x[i]` is aligned with a gap
/// ```ignore
///              .... A   G  x_i
///              .... G  y_j  -
/// ```
/// This is interpreted as an insertion into `x` w.r.t reference `y`
///
/// `D(i,j)` is the best score such that `y[j]` is aligned with a gap
/// ```ignore
///              .... A  x_i  -
///              .... G   G  y_j
/// ```
/// This is interpreted as a deletion from `x` w.r.t reference `y`
///
/// `S(i,j)` is the best score for prefixes `x[0..i]`, `y[0..j]`
///
/// To save space, only two columns of these matrices are stored at
/// any point - the current column and the previous one. Moreover
/// `M(i,j)` is not explicitly stored
///
/// `Lx` is the optimal x suffix clipping lengths from each position of the
/// sequence y
///
/// `Ly` is the optimal y suffix clipping lengths from each position of the
/// sequence x
///
/// `Sn` is the last column of the matrix. This is needed to keep track of
/// suffix clipping scores
///
/// `traceback` - see [`bio::alignment::pairwise::TracebackCell`](struct.TracebackCell.html)
///
/// `scoring` - see [`bio::alignment::pairwise::Scoring`](struct.Scoring.html)
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct DoubleStrandAligner<F: MatchFunc> {
    forward: SingleStrandMatrix<F>,
    reverse: SingleStrandMatrix<F>,
}

impl<F: MatchFunc> DoubleStrandAligner<F> {
    /// Create new aligner instance with given gap open and gap extend penalties
    /// and the score function.
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `jump_score` - the score for jumping back in the query (should not be positive)
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    pub fn new(gap_open: i32, gap_extend: i32, jump_score: i32, match_fn: F) -> Self {
        DoubleStrandAligner::with_capacity(
            DEFAULT_ALIGNER_CAPACITY,
            DEFAULT_ALIGNER_CAPACITY,
            gap_open,
            gap_extend,
            jump_score,
            match_fn,
        )
    }

    /// Create new aligner instance. The size hints help to
    /// avoid unnecessary memory allocations.
    ///
    /// # Arguments
    ///
    /// * `m` - the expected size of x
    /// * `n` - the expected size of y
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `jump_score` - the score for jumping back in the query (should not be positive)
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    pub fn with_capacity(
        m: usize,
        n: usize,
        gap_open: i32,
        gap_extend: i32,
        jump_score: i32,
        match_fn: F,
    ) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        DoubleStrandAligner {
            forward: SingleStrandMatrix::with_capacity(
                m, n, gap_open, gap_extend, jump_score, match_fn,
            ),
            reverse: SingleStrandMatrix::with_capacity(
                m, n, gap_open, gap_extend, jump_score, match_fn,
            ),
        }
    }

    /// Create new aligner instance with given the scoring struct
    ///
    /// # Arguments
    ///
    /// * `scoring` - the scoring struct (see bio::alignment::pairwise::Scoring)
    pub fn with_scoring(scoring: Scoring<F>) -> Self {
        DoubleStrandAligner::with_capacity_and_scoring(
            DEFAULT_ALIGNER_CAPACITY,
            DEFAULT_ALIGNER_CAPACITY,
            scoring,
        )
    }

    /// Create new aligner instance with scoring and size hint. The size hints help to
    /// avoid unnecessary memory allocations.
    ///
    /// # Arguments
    ///
    /// * `m` - the expected size of x
    /// * `n` - the expected size of y
    /// * `scoring` - the scoring struct
    pub fn with_capacity_and_scoring(m: usize, n: usize, scoring: Scoring<F>) -> Self {
        assert!(scoring.gap_open <= 0, "gap_open can't be positive");
        assert!(scoring.gap_extend <= 0, "gap_extend can't be positive");
        assert!(
            scoring.xclip_prefix <= 0,
            "Clipping penalty (x prefix) can't be positive"
        );
        assert!(
            scoring.xclip_suffix <= 0,
            "Clipping penalty (x suffix) can't be positive"
        );
        assert!(
            scoring.yclip_prefix <= 0,
            "Clipping penalty (y prefix) can't be positive"
        );
        assert!(
            scoring.yclip_suffix <= 0,
            "Clipping penalty (y suffix) can't be positive"
        );

        DoubleStrandAligner {
            forward: SingleStrandMatrix::with_capacity_and_scoring(m, n, scoring),
            reverse: SingleStrandMatrix::with_capacity_and_scoring(m, n, scoring),
        }
    }

    fn do_traceback(
        &mut self,
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        m: usize,
        n: usize,
    ) -> PairwiseAlignment {
        let mut i = m;
        let mut j = n;
        let mut operations: Vec<AlignmentOperation> = Vec::with_capacity(x.len());
        let mut xstart: usize = 0usize;
        let mut ystart: usize = 0usize;
        let mut xend = m;
        let mut yend = n;

        let mut on_forward: bool = self.forward.S[i][j] > self.reverse.S[i][j];
        let aligner = if on_forward {
            &self.forward
        } else {
            &self.reverse
        };

        let mut last_layer = aligner.traceback.get(i, j).get_s_bits();
        loop {
            let next_layer: u32;
            match last_layer {
                TB_START => break,
                TB_INS => {
                    operations.push(AlignmentOperation::Ins);
                    next_layer = aligner.traceback.get(i, j).get_i_bits();
                    i -= 1;
                }
                TB_DEL => {
                    operations.push(AlignmentOperation::Del);
                    next_layer = aligner.traceback.get(i, j).get_d_bits();
                    j -= 1;
                }
                TB_MATCH | TB_SUBST => {
                    if last_layer == TB_MATCH {
                        operations.push(AlignmentOperation::Match);
                    } else {
                        operations.push(AlignmentOperation::Subst);
                    }
                    let from_i = aligner.traceback.get(i, j).get_l_bits() as usize;
                    if from_i != i - 1 {
                        operations.push(AlignmentOperation::Xskip(i - 1));
                    }
                    next_layer = aligner.traceback.get(from_i, j - 1).get_s_bits();
                    i = from_i;
                    j -= 1;
                }
                TB_XCLIP_PREFIX => {
                    operations.push(AlignmentOperation::Xclip(i));
                    xstart = i;
                    i = 0;
                    next_layer = aligner.traceback.get(0, j).get_s_bits();
                }
                TB_XCLIP_SUFFIX => {
                    operations.push(AlignmentOperation::Xclip(aligner.Lx[j]));
                    i -= aligner.Lx[j];
                    xend = i;
                    next_layer = aligner.traceback.get(i, j).get_s_bits();
                }
                TB_YCLIP_PREFIX => {
                    operations.push(AlignmentOperation::Yclip(j));
                    ystart = j;
                    j = 0;
                    next_layer = aligner.traceback.get(i, 0).get_s_bits();
                }
                TB_YCLIP_SUFFIX => {
                    operations.push(AlignmentOperation::Yclip(aligner.Ly[i]));
                    j -= aligner.Ly[i];
                    yend = j;
                    next_layer = aligner.traceback.get(i, j).get_s_bits();
                }
                _ => panic!("Dint expect this!"),
            }
            last_layer = next_layer;
        }

        operations.reverse();
        PairwiseAlignment {
            score: aligner.S[n % 2][m],
            ystart,
            xstart,
            yend,
            xend,
            ylen: n,
            xlen: m,
            operations,
            mode: AlignmentMode::Custom,
        }
    }

    /// The core function to compute the alignment
    ///
    /// # Arguments
    ///
    /// * `x` - Textslice
    /// * `y` - Textslice
    pub fn custom(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> PairwiseAlignment {
        let (m, n) = (x.len(), y.len());

        // Set the initial conditions
        // We are repeating some work, but that's okay!
        self.forward.init_matrices(m, n);
        self.reverse.init_matrices(m, n);

        for j in 1..=n {
            let curr = j % 2;
            let prev = 1 - curr;

            // Initialize the column
            self.forward.init_column(j, prev, curr, m, n);
            self.reverse.init_column(j, prev, curr, m, n);

            // Fill the column
            self.forward.fill_column(x, y, m, n, j, prev, curr);
            self.reverse.fill_column(x, y, m, n, j, prev, curr);
        }

        self.forward.fill_last_column_and_end_clipping(m, n);
        self.reverse.fill_last_column_and_end_clipping(m, n);

        // Traceback...
        self.do_traceback(x, y, m, n)
    }

    /// Calculate global alignment of x against y.
    pub fn global(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> PairwiseAlignment {
        // Store the current clip penalties
        let forward_clip_penalties = [
            self.forward.scoring.xclip_prefix,
            self.forward.scoring.xclip_suffix,
            self.forward.scoring.yclip_prefix,
            self.forward.scoring.yclip_suffix,
        ];
        let reverse_clip_penalties = [
            self.reverse.scoring.xclip_prefix,
            self.reverse.scoring.xclip_suffix,
            self.reverse.scoring.yclip_prefix,
            self.reverse.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.forward.scoring.xclip_prefix = MIN_SCORE;
        self.forward.scoring.xclip_suffix = MIN_SCORE;
        self.forward.scoring.yclip_prefix = MIN_SCORE;
        self.forward.scoring.yclip_suffix = MIN_SCORE;
        self.reverse.scoring.xclip_prefix = MIN_SCORE;
        self.reverse.scoring.xclip_suffix = MIN_SCORE;
        self.reverse.scoring.yclip_prefix = MIN_SCORE;
        self.reverse.scoring.yclip_suffix = MIN_SCORE;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::Global;

        // Set the clip penalties to the original values
        self.forward.scoring.xclip_prefix = forward_clip_penalties[0];
        self.forward.scoring.xclip_suffix = forward_clip_penalties[1];
        self.forward.scoring.yclip_prefix = forward_clip_penalties[2];
        self.forward.scoring.yclip_suffix = forward_clip_penalties[3];
        self.reverse.scoring.xclip_prefix = reverse_clip_penalties[0];
        self.reverse.scoring.xclip_suffix = reverse_clip_penalties[1];
        self.reverse.scoring.yclip_prefix = reverse_clip_penalties[2];
        self.reverse.scoring.yclip_suffix = reverse_clip_penalties[3];
        alignment
    }

    /// Calculate semiglobal alignment of x against y (x is global, y is local).
    pub fn semiglobal(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> PairwiseAlignment {
        // Store the current clip penalties
        let forward_clip_penalties = [
            self.forward.scoring.xclip_prefix,
            self.forward.scoring.xclip_suffix,
            self.forward.scoring.yclip_prefix,
            self.forward.scoring.yclip_suffix,
        ];
        let reverse_clip_penalties = [
            self.reverse.scoring.xclip_prefix,
            self.reverse.scoring.xclip_suffix,
            self.reverse.scoring.yclip_prefix,
            self.reverse.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.forward.scoring.xclip_prefix = MIN_SCORE;
        self.forward.scoring.xclip_suffix = MIN_SCORE;
        self.forward.scoring.yclip_prefix = 0;
        self.forward.scoring.yclip_suffix = 0;
        self.reverse.scoring.xclip_prefix = MIN_SCORE;
        self.reverse.scoring.xclip_suffix = MIN_SCORE;
        self.reverse.scoring.yclip_prefix = 0;
        self.reverse.scoring.yclip_suffix = 0;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::Semiglobal;

        // Filter out Xclip and Yclip from alignment.operations
        {
            use self::AlignmentOperation::{Del, Ins, Match, Subst, Xclip, Xskip};
            alignment.operations.retain(|x| {
                *x == Match
                    || *x == Subst
                    || *x == Ins
                    || *x == Del
                    || matches!(*x, Xclip(_))
                    || matches!(*x, Xskip(_))
            });
        }

        // Set the clip penalties to the original values

        self.forward.scoring.xclip_prefix = forward_clip_penalties[0];
        self.forward.scoring.xclip_suffix = forward_clip_penalties[1];
        self.forward.scoring.yclip_prefix = forward_clip_penalties[2];
        self.forward.scoring.yclip_suffix = forward_clip_penalties[3];
        self.reverse.scoring.xclip_prefix = reverse_clip_penalties[0];
        self.reverse.scoring.xclip_suffix = reverse_clip_penalties[1];
        self.reverse.scoring.yclip_prefix = reverse_clip_penalties[2];
        self.reverse.scoring.yclip_suffix = reverse_clip_penalties[3];

        alignment
    }

    /// Calculate semiglobal alignment of x against y (x is local, y is global).
    pub fn semiglobal2(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> PairwiseAlignment {
        // Store the current clip penalties
        let forward_clip_penalties = [
            self.forward.scoring.xclip_prefix,
            self.forward.scoring.xclip_suffix,
            self.forward.scoring.yclip_prefix,
            self.forward.scoring.yclip_suffix,
        ];
        let reverse_clip_penalties = [
            self.reverse.scoring.xclip_prefix,
            self.reverse.scoring.xclip_suffix,
            self.reverse.scoring.yclip_prefix,
            self.reverse.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.forward.scoring.xclip_prefix = 0;
        self.forward.scoring.xclip_suffix = 0;
        self.forward.scoring.yclip_prefix = MIN_SCORE;
        self.forward.scoring.yclip_suffix = MIN_SCORE;
        self.reverse.scoring.xclip_prefix = 0;
        self.reverse.scoring.xclip_suffix = 0;
        self.reverse.scoring.yclip_prefix = MIN_SCORE;
        self.reverse.scoring.yclip_suffix = MIN_SCORE;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::Semiglobal;

        // Filter out Xclip and Yclip from alignment.operations
        {
            use self::AlignmentOperation::{Del, Ins, Match, Subst, Xclip, Xskip};
            alignment.operations.retain(|x| {
                *x == Match
                    || *x == Subst
                    || *x == Ins
                    || *x == Del
                    || matches!(*x, Xclip(_))
                    || matches!(*x, Xskip(_))
            });
        }

        // Set the clip penalties to the original values

        self.forward.scoring.xclip_prefix = forward_clip_penalties[0];
        self.forward.scoring.xclip_suffix = forward_clip_penalties[1];
        self.forward.scoring.yclip_prefix = forward_clip_penalties[2];
        self.forward.scoring.yclip_suffix = forward_clip_penalties[3];
        self.reverse.scoring.xclip_prefix = reverse_clip_penalties[0];
        self.reverse.scoring.xclip_suffix = reverse_clip_penalties[1];
        self.reverse.scoring.yclip_prefix = reverse_clip_penalties[2];
        self.reverse.scoring.yclip_suffix = reverse_clip_penalties[3];

        alignment
    }

    /// Calculate local alignment of x against y.
    pub fn local(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> PairwiseAlignment {
        // Store the current clip penalties
        let forward_clip_penalties = [
            self.forward.scoring.xclip_prefix,
            self.forward.scoring.xclip_suffix,
            self.forward.scoring.yclip_prefix,
            self.forward.scoring.yclip_suffix,
        ];
        let reverse_clip_penalties = [
            self.reverse.scoring.xclip_prefix,
            self.reverse.scoring.xclip_suffix,
            self.reverse.scoring.yclip_prefix,
            self.reverse.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.forward.scoring.xclip_prefix = 0;
        self.forward.scoring.xclip_suffix = 0;
        self.forward.scoring.yclip_prefix = 0;
        self.forward.scoring.yclip_suffix = 0;
        self.reverse.scoring.xclip_prefix = 0;
        self.reverse.scoring.xclip_suffix = 0;
        self.reverse.scoring.yclip_prefix = 0;
        self.reverse.scoring.yclip_suffix = 0;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::Local;

        // Filter out Xclip and Yclip from alignment.operations
        {
            use self::AlignmentOperation::{Del, Ins, Match, Subst, Xclip, Xskip};
            alignment.operations.retain(|x| {
                *x == Match
                    || *x == Subst
                    || *x == Ins
                    || *x == Del
                    || matches!(*x, Xclip(_))
                    || matches!(*x, Xskip(_))
            });
        }

        // Set the clip penalties to the original values

        self.forward.scoring.xclip_prefix = forward_clip_penalties[0];
        self.forward.scoring.xclip_suffix = forward_clip_penalties[1];
        self.forward.scoring.yclip_prefix = forward_clip_penalties[2];
        self.forward.scoring.yclip_suffix = forward_clip_penalties[3];
        self.reverse.scoring.xclip_prefix = reverse_clip_penalties[0];
        self.reverse.scoring.xclip_suffix = reverse_clip_penalties[1];
        self.reverse.scoring.yclip_prefix = reverse_clip_penalties[2];
        self.reverse.scoring.yclip_suffix = reverse_clip_penalties[3];

        alignment
    }
}
