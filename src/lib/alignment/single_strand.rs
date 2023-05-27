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
use crate::alignment::scoring::Scoring;
use bio::alignment::pairwise::MatchFunc;
use bio::alignment::pairwise::MatchParams;
use bio::utils::TextSlice;

use crate::alignment::constants::AlignmentOperation;
use crate::alignment::constants::DEFAULT_ALIGNER_CAPACITY;
use crate::alignment::pairwise::PairwiseAlignment;
use crate::alignment::traceback::Traceback;
use crate::alignment::traceback::TracebackCell;
use crate::alignment::x_buffer::XBuffer;

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
/// (or backward) in `x` (on the same strand).
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
#[allow(non_snake_case)]
pub struct SingleStrandAligner<F: MatchFunc> {
    pub I: [Vec<i32>; 2],
    pub D: [Vec<i32>; 2],
    pub S: [Vec<i32>; 2],
    pub Lx: Vec<usize>,
    pub Ly: Vec<usize>,
    pub Sn: Vec<i32>,
    pub traceback: Traceback,
    pub scoring: Scoring<F>,
    pub x_buffer: XBuffer,
}

impl Default for SingleStrandAligner<MatchParams> {
    fn default() -> Self {
        let match_fn = MatchParams::new(1, -1);
        SingleStrandAligner::new(-5, -1, -10, match_fn)
    }
}

impl<F: MatchFunc> SingleStrandAligner<F> {
    pub fn init_matrices(&mut self, m: usize, n: usize) {
        // initialize the traceback
        self.traceback.init(m, n);

        // Set the initial conditions
        // We are repeating some work, but that's okay!
        for k in 0..2 {
            self.I[k].clear();
            self.D[k].clear();
            self.S[k].clear();

            self.D[k].extend(repeat(MIN_SCORE).take(m + 1));
            self.I[k].extend(repeat(MIN_SCORE).take(m + 1));
            self.S[k].extend(repeat(MIN_SCORE).take(m + 1));

            self.S[k][0] = 0;

            if k == 0 {
                let mut tb = TracebackCell::new();
                tb.set_all(TB_START);
                tb.set_s_flip_strand(false);
                tb.set_i_flip_strand(false);
                tb.set_s_from(0);
                tb.set_i_from(0);
                self.traceback.set(0, 0, tb);
                self.Lx.clear();
                self.Lx.extend(repeat(0usize).take(n + 1));
                self.Ly.clear();
                self.Ly.extend(repeat(0usize).take(m + 1));
                self.Sn.clear();
                self.Sn.extend(repeat(MIN_SCORE).take(m + 1));
                self.Sn[0] = self.scoring.yclip_suffix;
                self.Ly[0] = n;
            }

            for i in 1..=m {
                let mut tb = TracebackCell::new();
                tb.set_all(TB_START);
                tb.set_s_flip_strand(false);
                tb.set_i_flip_strand(false);
                tb.set_s_from(0);
                tb.set_i_from(0);
                if i == 1 {
                    self.I[k][i] = self.scoring.gap_open + self.scoring.gap_extend;
                    tb.set_i_bits(TB_START);
                    tb.set_i_from((i - 1) as u32);
                } else {
                    // Insert all i characters
                    let i_score = self.scoring.gap_open + self.scoring.gap_extend * (i as i32);
                    let c_score =
                        self.scoring.xclip_prefix + self.scoring.gap_open + self.scoring.gap_extend; // Clip then insert
                    if i_score > c_score {
                        self.I[k][i] = i_score;
                        tb.set_i_bits(TB_INS);
                        tb.set_i_from((i - 1) as u32);
                    } else {
                        self.I[k][i] = c_score;
                        tb.set_i_bits(TB_XCLIP_PREFIX);
                        tb.set_i_from((i - 1) as u32);
                    }
                }

                if i == m {
                    tb.set_s_bits(TB_XCLIP_SUFFIX);
                } else {
                    self.S[k][i] = MIN_SCORE;
                }

                if self.I[k][i] > self.S[k][i] {
                    self.S[k][i] = self.I[k][i];
                    tb.set_s_bits(TB_INS);
                }

                if self.scoring.xclip_prefix > self.S[k][i] {
                    self.S[k][i] = self.scoring.xclip_prefix;
                    tb.set_s_bits(TB_XCLIP_PREFIX);
                }

                // Track the score if we do a suffix clip (x) after this character
                if i != m && self.S[k][i] + self.scoring.xclip_suffix > self.S[k][m] {
                    self.S[k][m] = self.S[k][i] + self.scoring.xclip_suffix;
                    self.Lx[0] = m - i;
                }

                if k == 0 {
                    self.traceback.set(i, 0, tb);
                }
                // Track the score if we do suffix clip (y) from here
                if self.S[k][i] + self.scoring.yclip_suffix > self.Sn[i] {
                    self.Sn[i] = self.S[k][i] + self.scoring.yclip_suffix;
                    self.Ly[i] = n;
                }
            }
        }
    }

    pub fn init_column(&mut self, j: usize, curr: usize, m: usize, n: usize) {
        // Handle i = 0 case
        let mut tb = TracebackCell::new();
        tb.set_s_flip_strand(false);
        tb.set_i_flip_strand(false);
        tb.set_s_from(0);
        tb.set_i_from(0);
        self.I[curr][0] = MIN_SCORE;

        if j == 1 {
            self.D[curr][0] = self.scoring.gap_open + self.scoring.gap_extend;
            tb.set_d_bits(TB_START);
        } else {
            // Delete all j characters
            let d_score = self.scoring.gap_open + self.scoring.gap_extend * (j as i32);
            let c_score =
                self.scoring.yclip_prefix + self.scoring.gap_open + self.scoring.gap_extend;
            if d_score > c_score {
                self.D[curr][0] = d_score;
                tb.set_d_bits(TB_DEL);
            } else {
                self.D[curr][0] = c_score;
                tb.set_d_bits(TB_YCLIP_PREFIX);
            }
        }
        if self.D[curr][0] > self.scoring.yclip_prefix {
            self.S[curr][0] = self.D[curr][0];
            tb.set_s_bits(TB_DEL);
        } else {
            self.S[curr][0] = self.scoring.yclip_prefix;
            tb.set_s_bits(TB_YCLIP_PREFIX);
        }

        if j == n && self.Sn[0] > self.S[curr][0] {
            // Check if the suffix clip score is better
            self.S[curr][0] = self.Sn[0];
            tb.set_s_bits(TB_YCLIP_SUFFIX);
        // Track the score if we do suffix clip (y) from here
        } else if self.S[curr][0] + self.scoring.yclip_suffix > self.Sn[0] {
            self.Sn[0] = self.S[curr][0] + self.scoring.yclip_suffix;
            self.Ly[0] = n - j;
        }

        self.traceback.set(0, j, tb);

        for i in 1..=m {
            self.S[curr][i] = MIN_SCORE;
        }
    }

    pub fn fill_column(
        &mut self,
        x: TextSlice<'_>,
        y: TextSlice<'_>,
        m: usize,
        n: usize,
        j: usize,
        prev: usize,
        curr: usize,
    ) {
        let q = y[j - 1];
        let xclip_score = self.scoring.xclip_prefix
            + max(
                self.scoring.yclip_prefix,
                self.scoring.gap_open + self.scoring.gap_extend * (j as i32),
            );
        for i in 1..=m {
            let p: u8 = x[i - 1];
            let mut tb = TracebackCell::new();

            // Align the x[i-1] with y[j-1], either through a jump or a diagonal move.
            let (x_jump_score, x_jump_from, x_jump_flip_strand) = self.x_buffer.get(i);
            let diag_score = self.S[prev][i - 1];
            let (m_score, m_from, m_flip_strand) = {
                let addend = self.scoring.match_fn.score(p, q);
                if diag_score >= x_jump_score {
                    (diag_score + addend, (i - 1) as u32, false)
                } else {
                    (x_jump_score + addend, x_jump_from, x_jump_flip_strand)
                }
            };

            // Align an insertion, either through moving from the previous base in x (gap open or
            // extend), or through a jump (gap open).
            let i_score = self.I[curr][i - 1] + self.scoring.gap_extend;
            let s_score = self.S[curr][i - 1] + self.scoring.gap_open + self.scoring.gap_extend;
            let j_score = x_jump_score + self.scoring.gap_open + self.scoring.gap_extend;
            let best_i_score = max(i_score, max(s_score, j_score));
            if i_score == best_i_score {
                tb.set_i_bits(TB_INS);
                tb.set_i_flip_strand(false);
                tb.set_i_from((i - 1) as u32);
            } else if s_score == best_i_score {
                tb.set_i_bits(self.traceback.get(i - 1, j).get_s_bits());
                tb.set_i_flip_strand(false);
                tb.set_i_from((i - 1) as u32);
            } else {
                tb.set_i_bits(self.traceback.get(x_jump_from as usize, j).get_s_bits());
                tb.set_i_flip_strand(x_jump_flip_strand);
                tb.set_i_from(x_jump_from);
            }

            let d_score = self.D[prev][i] + self.scoring.gap_extend;
            let s_score = self.S[prev][i] + self.scoring.gap_open + self.scoring.gap_extend;
            let best_d_score;
            if d_score > s_score {
                best_d_score = d_score;
                tb.set_d_bits(TB_DEL);
            } else {
                best_d_score = s_score;
                tb.set_d_bits(self.traceback.get(i, j - 1).get_s_bits());
            }

            tb.set_s_bits(TB_XCLIP_SUFFIX);
            let mut best_s_score = self.S[curr][i];

            if m_score > best_s_score {
                best_s_score = m_score;
                tb.set_s_bits(if p == q { TB_MATCH } else { TB_SUBST });
                tb.set_s_flip_strand(m_flip_strand);
                tb.set_s_from(m_from);
            }

            if best_d_score > best_s_score {
                best_s_score = best_d_score;
                tb.set_s_bits(TB_DEL);
            }

            if best_i_score > best_s_score {
                best_s_score = best_i_score;
                tb.set_s_bits(TB_INS);
            }

            if xclip_score > best_s_score {
                best_s_score = xclip_score;
                tb.set_s_bits(TB_XCLIP_PREFIX);
            }

            let yclip_score = self.scoring.yclip_prefix
                + self.scoring.gap_open
                + self.scoring.gap_extend * (i as i32);
            if yclip_score > best_s_score {
                best_s_score = yclip_score;
                tb.set_s_bits(TB_YCLIP_PREFIX);
            }

            self.S[curr][i] = best_s_score;
            self.I[curr][i] = best_i_score;
            self.D[curr][i] = best_d_score;

            // Track the score if we do suffix clip (x) from here
            if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][m] {
                self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
                self.Lx[j] = m - i;
            }

            // Track the score if we do suffix clip (y) from here
            if self.S[curr][i] + self.scoring.yclip_suffix > self.Sn[i] {
                self.Sn[i] = self.S[curr][i] + self.scoring.yclip_suffix;
                self.Ly[i] = n - j;
            }

            self.traceback.set(i, j, tb);
        }
    }

    pub fn fill_last_column_and_end_clipping(&mut self, m: usize, n: usize) {
        // Handle suffix clipping in the j=n case
        for i in 0..=m {
            let j = n;
            let curr = j % 2;
            if self.Sn[i] > self.S[curr][i] {
                self.S[curr][i] = self.Sn[i];
                self.traceback.get_mut(i, j).set_s_bits(TB_YCLIP_SUFFIX);
            }
            if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][m] {
                self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
                self.Lx[j] = m - i;
                self.traceback.get_mut(m, j).set_s_bits(TB_XCLIP_SUFFIX);
            }
        }

        // Since there could be a change in the last column of S,
        // recompute the last column of I as this could also change
        for i in 1..=m {
            let j = n;
            let curr = j % 2;
            let s_score = self.S[curr][i - 1] + self.scoring.gap_open + self.scoring.gap_extend;
            if s_score > self.I[curr][i] {
                self.I[curr][i] = s_score;
                let s_bit = self.traceback.get(i - 1, j).get_s_bits();
                self.traceback.get_mut(i, j).set_i_bits(s_bit);
            }
            if s_score > self.S[curr][i] {
                self.S[curr][i] = s_score;
                self.traceback.get_mut(i, j).set_s_bits(TB_INS);
                if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][m] {
                    self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
                    self.Lx[j] = m - i;
                    self.traceback.get_mut(m, j).set_s_bits(TB_XCLIP_SUFFIX);
                }
            }
        }
    }

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
        SingleStrandAligner::with_capacity(
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

        SingleStrandAligner {
            I: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            D: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            S: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            Lx: Vec::with_capacity(n + 1),
            Ly: Vec::with_capacity(m + 1),
            Sn: Vec::with_capacity(m + 1),
            traceback: Traceback::with_capacity(m, n),
            scoring: Scoring::new(gap_open, gap_extend, jump_score, match_fn),
            x_buffer: XBuffer::new(m),
        }
    }

    /// Create new aligner instance with given the scoring struct
    ///
    /// # Arguments
    ///
    /// * `scoring` - the scoring struct (see bio::alignment::pairwise::Scoring)
    pub fn with_scoring(scoring: Scoring<F>) -> Self {
        SingleStrandAligner::with_capacity_and_scoring(
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

        SingleStrandAligner {
            I: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            D: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            S: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            Lx: Vec::with_capacity(n + 1),
            Ly: Vec::with_capacity(m + 1),
            Sn: Vec::with_capacity(m + 1),
            traceback: Traceback::with_capacity(m, n),
            scoring,
            x_buffer: XBuffer::new(m),
        }
    }

    fn do_traceback(&mut self, m: usize, n: usize) -> PairwiseAlignment {
        let mut i = m;
        let mut j = n;
        let mut operations: Vec<AlignmentOperation> = Vec::with_capacity(m);
        let mut xstart: usize = 0usize;
        let mut ystart: usize = 0usize;
        let mut xend = m;
        let mut yend = n;

        let mut last_layer = self.traceback.get(i, j).get_s_bits();
        loop {
            let next_layer: u16;
            match last_layer {
                TB_START => break,
                TB_INS => {
                    operations.push(AlignmentOperation::Ins);
                    let i_from = self.traceback.get(i, j).get_i_from() as usize;
                    let i_flip_strand = self.traceback.get(i, j).get_i_flip_strand();
                    if i_from == i - 1 {
                        next_layer = self.traceback.get(i, j).get_i_bits();
                    } else {
                        if i_flip_strand {
                            operations.push(AlignmentOperation::Xflip(i_from));
                        } else {
                            operations.push(AlignmentOperation::Xskip(i_from));
                        }
                        next_layer = self.traceback.get(i_from, j).get_i_bits();
                    }
                    i = i_from;
                }
                TB_DEL => {
                    operations.push(AlignmentOperation::Del);
                    next_layer = self.traceback.get(i, j).get_d_bits();
                    j -= 1;
                }
                TB_MATCH | TB_SUBST => {
                    if last_layer == TB_MATCH {
                        operations.push(AlignmentOperation::Match);
                    } else {
                        operations.push(AlignmentOperation::Subst);
                    }
                    let s_from = self.traceback.get(i, j).get_s_from() as usize;
                    let s_flip_strand = self.traceback.get(i, j).get_s_flip_strand();
                    assert!(!s_flip_strand);
                    if s_from == i - 1 {
                        next_layer = self.traceback.get(i - 1, j - 1).get_s_bits();
                    } else {
                        operations.push(AlignmentOperation::Xskip(s_from));
                        next_layer = self.traceback.get(s_from, j - 1).get_s_bits();
                    }
                    i = s_from;
                    j -= 1;
                }
                TB_XCLIP_PREFIX => {
                    operations.push(AlignmentOperation::Xclip(i));
                    xstart = i;
                    i = 0;
                    next_layer = self.traceback.get(0, j).get_s_bits();
                }
                TB_XCLIP_SUFFIX => {
                    operations.push(AlignmentOperation::Xclip(self.Lx[j]));
                    i -= self.Lx[j];
                    xend = i;
                    next_layer = self.traceback.get(i, j).get_s_bits();
                }
                TB_YCLIP_PREFIX => {
                    operations.push(AlignmentOperation::Yclip(j));
                    ystart = j;
                    j = 0;
                    next_layer = self.traceback.get(i, 0).get_s_bits();
                }
                TB_YCLIP_SUFFIX => {
                    operations.push(AlignmentOperation::Yclip(self.Ly[i]));
                    j -= self.Ly[i];
                    yend = j;
                    next_layer = self.traceback.get(i, j).get_s_bits();
                }
                _ => panic!("Dint expect this!"),
            }
            last_layer = next_layer;
        }

        operations.reverse();
        PairwiseAlignment {
            score: self.S[n % 2][m],
            ystart,
            xstart,
            yend,
            xend,
            ylen: n,
            xlen: m,
            is_forward: true,
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

        self.init_matrices(m, n);

        for j in 1..=n {
            let curr = j % 2;
            let prev = 1 - curr;

            // Initialize the column
            self.init_column(j, curr, m, n);

            // Initiliaze the jump buffers
            self.x_buffer.fill(m, prev, &self.S, &self.scoring);

            // Fill the column
            self.fill_column(x, y, m, n, j, prev, curr);
        }

        self.fill_last_column_and_end_clipping(m, n);

        self.do_traceback(m, n)
    }

    /// Calculate global alignment of x against y.
    pub fn global(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> PairwiseAlignment {
        // Store the current clip penalties
        let clip_penalties = [
            self.scoring.xclip_prefix,
            self.scoring.xclip_suffix,
            self.scoring.yclip_prefix,
            self.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.scoring.xclip_prefix = MIN_SCORE;
        self.scoring.xclip_suffix = MIN_SCORE;
        self.scoring.yclip_prefix = MIN_SCORE;
        self.scoring.yclip_suffix = MIN_SCORE;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::Global;

        // Set the clip penalties to the original values
        self.scoring.xclip_prefix = clip_penalties[0];
        self.scoring.xclip_suffix = clip_penalties[1];
        self.scoring.yclip_prefix = clip_penalties[2];
        self.scoring.yclip_suffix = clip_penalties[3];

        alignment
    }

    /// Calculate semiglobal alignment of x against y (x is global, y is local).
    pub fn semiglobal(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> PairwiseAlignment {
        // Store the current clip penalties
        let clip_penalties = [
            self.scoring.xclip_prefix,
            self.scoring.xclip_suffix,
            self.scoring.yclip_prefix,
            self.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.scoring.xclip_prefix = MIN_SCORE;
        self.scoring.xclip_suffix = MIN_SCORE;
        self.scoring.yclip_prefix = 0;
        self.scoring.yclip_suffix = 0;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::Semiglobal;

        // Filter out Yclip from alignment.operations
        {
            use self::AlignmentOperation::{Del, Ins, Match, Subst, Xskip};
            alignment.operations.retain(|x| {
                *x == Match || *x == Subst || *x == Ins || *x == Del || matches!(*x, Xskip(_))
            });
        }

        // Set the clip penalties to the original values
        self.scoring.xclip_prefix = clip_penalties[0];
        self.scoring.xclip_suffix = clip_penalties[1];
        self.scoring.yclip_prefix = clip_penalties[2];
        self.scoring.yclip_suffix = clip_penalties[3];

        alignment
    }

    /// Calculate local alignment of x against y.
    pub fn local(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> PairwiseAlignment {
        // Store the current clip penalties
        let clip_penalties = [
            self.scoring.xclip_prefix,
            self.scoring.xclip_suffix,
            self.scoring.yclip_prefix,
            self.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.scoring.xclip_prefix = 0;
        self.scoring.xclip_suffix = 0;
        self.scoring.yclip_prefix = 0;
        self.scoring.yclip_suffix = 0;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::Local;

        // Filter out Xclip and Yclip from alignment.operations
        {
            use self::AlignmentOperation::{Del, Ins, Match, Subst, Xskip};
            alignment.operations.retain(|x| {
                *x == Match || *x == Subst || *x == Ins || *x == Del || matches!(*x, Xskip(_))
            });
        }

        // Set the clip penalties to the original values
        self.scoring.xclip_prefix = clip_penalties[0];
        self.scoring.xclip_suffix = clip_penalties[1];
        self.scoring.yclip_prefix = clip_penalties[2];
        self.scoring.yclip_suffix = clip_penalties[3];

        alignment
    }
}

// Tests
#[cfg(test)]
pub mod tests {
    use bio::alignment::pairwise::{MatchFunc, MatchParams};
    use itertools::Itertools;
    use rstest::rstest;

    use crate::alignment::constants::AlignmentOperation::{Del, Ins, Match, Subst, Xclip, Xskip};
    use crate::alignment::{constants::AlignmentOperation, pairwise::PairwiseAlignment};

    use super::SingleStrandAligner;

    /// Upper-cases and remove display-related characters from a string.
    fn s(bases: &str) -> Vec<u8> {
        bases
            .chars()
            .filter(|base| *base != '-' && *base != ' ' && *base != '_')
            .map(|base| base.to_ascii_uppercase() as u8)
            .collect_vec()
    }

    fn assert_alignment(
        alignment: &PairwiseAlignment,
        xstart: usize,
        xend: usize,
        ystart: usize,
        yend: usize,
        score: i32,
        operations: &[AlignmentOperation],
    ) {
        assert_eq!(alignment.xstart, xstart, "xstart {alignment:?}");
        assert_eq!(alignment.xend, xend, "xend {alignment:?}");
        assert_eq!(alignment.ystart, ystart, "ystart {alignment:?}");
        assert_eq!(alignment.yend, yend, "yend {alignment:?}");
        assert_eq!(alignment.score, score, "score {alignment:?}");
        assert_eq!(alignment.operations, operations, "operations {alignment:?}");
    }

    /// Identical sequences, all matches
    #[rstest]
    fn test_identical() {
        let x = s("ACGTAACC");
        let y = s("ACGTAACC");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8, &[Match; 8]);
    }

    /// Identical sequences, except one mismatch
    #[rstest]
    fn test_single_mismatch() {
        let x = s("AACCGGTT");
        let y = s("AACCGtTT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            8,
            0,
            8,
            7 - 1,
            &[Match, Match, Match, Match, Match, Subst, Match, Match],
        );
    }

    /// Identical sequences, except one samll deletion
    #[rstest]
    fn test_small_deletion() {
        let x = s("AACC-GTT");
        let y = s("AACCGGTT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            7,
            0,
            8,
            7 - (5 + 1),
            &[Match, Match, Match, Match, Del, Match, Match, Match],
        );
    }

    /// Identical sequences, except one samll insertion
    #[rstest]
    fn test_small_insertion() {
        let x = s("AACCGGTT");
        let y = s("AACC-GTT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            8,
            0,
            7,
            7 - (5 + 1),
            &[Match, Match, Match, Match, Ins, Match, Match, Match],
        );
    }

    /// Two sequences with compensating insertions and deletions
    #[rstest]
    fn test_compensating_insertion_and_deletion() {
        let x = s("AAACGCGCGCGCG-TT");
        let y = s("-AACGCGCGCGCGTTT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            15,
            0,
            15,
            14 - (5 + 1) - (5 + 1),
            &[
                Ins, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match,
                Match, Del, Match, Match,
            ],
        );
    }

    #[rstest]
    fn test_leading_insertion() {
        let x = s("ATTTTTTTTTTT");
        let y = s("-TTTTTTTTTTT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            12,
            0,
            11,
            11 - (5 + 1),
            &[
                Ins, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match,
            ],
        );
    }

    #[rstest]
    fn test_trailing_insertion() {
        let x = s("TTTTTTTTTTTA");
        let y = s("TTTTTTTTTTT-");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            12,
            0,
            11,
            11 - (5 + 1),
            &[
                Match, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match, Ins,
            ],
        );
    }

    #[rstest]
    fn test_leading_deletion() {
        let x = s("-TTTTTTTTTTT");
        let y = s("ATTTTTTTTTTT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            11,
            0,
            12,
            11 - (5 + 1),
            &[
                Del, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match,
            ],
        );
    }

    #[rstest]
    fn test_trailing_deletion() {
        let x = s("TTTTTTTTTTT-");
        let y = s("TTTTTTTTTTTA");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            11,
            0,
            12,
            11 - (5 + 1),
            &[
                Match, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match, Del,
            ],
        );
    }

    /// Align two sequences preferring a 2bp insertion and mismatch vs two small insertions"
    #[rstest]
    fn test_prefer_2bp_insertion_and_mismatch_vs_two_small_insertions() {
        let x = s("ATTTTTTTTTTTA");
        let y = s("--TTTTTTTTTTt");
        let match_fn = MatchParams::new(1, -1);
        let mut aligner = SingleStrandAligner::new(-3, -1, -10, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            13,
            0,
            11,
            10 - (3 + 1) - 1 - 1,
            &[
                Ins, Ins, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match,
                Subst,
            ],
        );
    }

    /// Align two sequences preferring a 2bp insertion and mismatch vs two small insertions"
    #[rstest]
    fn test_prefer_two_small_insertions_vs_2bp_insertion_and_mismatch() {
        let x = s("ATTTTTTTTTTTA");
        let y = s("-TTTTTTTTTTT-");
        let match_fn = MatchParams::new(1, -3);
        let mut aligner = SingleStrandAligner::new(-3, -1, -10, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            13,
            0,
            11,
            11 - (3 + 1) - (3 + 1),
            &[
                Ins, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match,
                Ins,
            ],
        );
    }

    #[rstest]
    fn test_left_justify_insertion_in_homopolymer() {
        let x = s("GTTTTTTTTTTA");
        let y = s("G-TTTTTTTTTA");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            12,
            0,
            11,
            11 - (5 + 1),
            &[
                Match, Ins, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match,
            ],
        );
    }

    #[rstest]
    fn test_left_justify_insertion_in_triplet() {
        let x = s("GACGACGACGACGA");
        let y = s("---GACGACGACGA");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            14,
            0,
            11,
            11 - (5 + 1) - 1 - 1,
            &[
                Ins, Ins, Ins, Match, Match, Match, Match, Match, Match, Match, Match, Match,
                Match, Match,
            ],
        );
    }

    #[rstest]
    fn test_left_justify_insertion_in_triplet_with_leading_matches() {
        let x = s("TTTGACGACGACGACGA");
        let y = s("TTT---GACGACGACGA");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            17,
            0,
            14,
            14 - (5 + 1) - 1 - 1,
            &[
                Match, Match, Match, Ins, Ins, Ins, Match, Match, Match, Match, Match, Match,
                Match, Match, Match, Match, Match,
            ],
        );
    }

    /// NB: if the jump score is set to -10, then the alignment would jump back in x (3bp), to give
    /// 17 matches (+17) and one jump (-10) for a score of 7, which is greater than 14 matches (+14)
    /// and a 3bp deletion (-5 - 1 - 1 -1 = -8) for a score of 6.
    #[rstest]
    fn test_left_justify_deletion_in_triplet_with_leading_matches() {
        let x = s("TTT---GACGACGACGA");
        let y = s("TTTGACGACGACGACGA");
        let mut aligner = SingleStrandAligner::default();
        aligner.scoring.jump_score = -11;
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            14,
            0,
            17,
            14 - (5 + 1) - 1 - 1,
            &[
                Match, Match, Match, Del, Del, Del, Match, Match, Match, Match, Match, Match,
                Match, Match, Match, Match, Match,
            ],
        );
    }

    /// NB: if the jump score is set to -10, then the alignment would jump back in x (3bp), to give
    /// 17 matches (+17) and one jump (-10) for a score of 7, which is greater than 14 matches (+14)
    /// and a 3bp deletion (-5 - 1 - 1 -1 = -8) for a score of 6.
    #[rstest]
    fn test_x_skip_in_triplet_with_leading_matches() {
        let x = s("TTT---GACGACGACGA");
        let y = s("TTTGACGACGACGACGA");
        let mut aligner = SingleStrandAligner::default();
        aligner.scoring.jump_score = -10;
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            14,
            0,
            17,
            17 - 10,
            &[
                Match,
                Match,
                Match,
                Match,
                Match,
                Match,
                Xskip(6),
                Match,
                Match,
                Match,
                Match,
                Match,
                Match,
                Match,
                Match,
                Match,
                Match,
                Match,
            ],
        );
    }

    /// Prefer a mismatch over an insertion and deletion when the mismatch score is less than than
    /// two times the gap open + gap extend scores. NB: could be either "1I2=1D3=" or "2=1X3=".
    #[rstest]
    fn test_prefer_mismatch_over_insertion_and_deletion() {
        let x = s("AAACCC");
        let y = s("AAcCCC");
        let match_fn = MatchParams::new(1, -3);
        let mut aligner = SingleStrandAligner::new(-1, -1, -10, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            6,
            0,
            6,
            5 - 3,
            &[Match, Match, Subst, Match, Match, Match],
        );
    }

    #[rstest]
    fn test_prefer_mismatch_over_insertion_and_deletion_when_same_score() {
        let x = s("AAACCC");
        let y = s("AAcCCC");
        // NB: could be either "1I2=1D3=" or "2=1X3="
        let match_fn = MatchParams::new(1, -4);
        let mut aligner = SingleStrandAligner::new(-1, -1, -10, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            6,
            0,
            6,
            5 - 4,
            &[Match, Match, Subst, Match, Match, Match],
        );
    }

    /// Prefer an insertion and deletion over a mismatch when the mismatch score is greater than
    /// two times the gap open + gap extend scores.
    #[rstest]
    fn test_prefer_insertion_and_deletion_over_mismatch() {
        let x = s("AAA-CCC");
        let y = s("-AACCCC");
        // NB: could be either "1I2=1D3=" or "2=1X3="
        // NB: if we prefer a insertion over a deletion, then it would be 2M1D1I3M
        let match_fn = MatchParams::new(1, -5);
        let mut aligner = SingleStrandAligner::new(-1, -1, -10, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            6,
            0,
            6,
            5 - (1 + 1) - (1 + 1),
            &[Ins, Match, Match, Del, Match, Match, Match],
        );
    }

    #[rstest]
    fn test_prefer_one_insertion_with_large_gap_open() {
        let x = s("ATTTTTTTTTTTA");
        let y = s("--TTTTTTTTTTt");
        let match_fn = MatchParams::new(1, -5);
        let mut aligner = SingleStrandAligner::new(-100, -1, -10000, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            13,
            0,
            11,
            10 - (100 + 1) - 1 - 5,
            &[
                Ins, Ins, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match,
                Subst,
            ],
        );
    }

    #[rstest]
    fn test_prefer_two_small_insertions_with_large_gap_extend() {
        let x = s("ATTTTTTTTTTTA");
        let y = s("-TTTTTTTTTTT-");
        let match_fn = MatchParams::new(1, -5);
        let mut aligner = SingleStrandAligner::new(-1, -100, -10000, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            13,
            0,
            11,
            11 - (100 + 1) - (100 + 1),
            &[
                Ins, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match,
                Ins,
            ],
        );
    }

    #[rstest]
    fn test_semiglobal_identical() {
        let x = s("ACGTAACC");
        let y = s("ACGTAACC");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.semiglobal(&x, &y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8, &[Match; 8]);
    }

    #[rstest]
    fn test_semiglobal_identical_subsequence() {
        let x = s("  CCGG  ");
        let y = s("AACCGGTT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.semiglobal(&x, &y);
        assert_alignment(&alignment, 0, 4, 2, 6, 4, &[Match; 4]);
    }

    #[rstest]
    fn test_semiglobal_subsequence_with_mismatch() {
        let x = s("       CGCGTCGTATACGTCGTT");
        let y = s("AAGATATCGCGTCGTATACGTCGTa");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.semiglobal(&x, &y);
        assert_alignment(
            &alignment,
            0,
            18,
            7,
            25,
            17 - 1,
            &[
                Match, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match, Match,
                Match, Match, Match, Match, Match, Subst,
            ],
        );
    }

    #[rstest]
    fn test_semiglobal_subsequence_with_deletion() {
        let x = s("  CGCG-CGCG  ");
        let y = s("AACGCGACGCGTT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.semiglobal(&x, &y);
        assert_alignment(
            &alignment,
            0,
            8,
            2,
            11,
            8 - (5 + 1),
            &[Match, Match, Match, Match, Del, Match, Match, Match, Match],
        );
    }

    #[rstest]
    fn test_semiglobal_insertion_when_x_longer_than_y() {
        let x = s("AAAAGGGGTTTT");
        let y = s("AAAA----TTTT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.semiglobal(&x, &y);
        assert_alignment(
            &alignment,
            0,
            12,
            0,
            8,
            8 - (5 + 1) - 1 - 1 - 1,
            &[
                Match, Match, Match, Match, Ins, Ins, Ins, Ins, Match, Match, Match, Match,
            ],
        );
    }

    #[rstest]
    fn test_global_leading_and_trailing_deletions() {
        let x = s("-------------------GGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTG---------------------------");
        let y = s("AGGGCTATAGACTGCTAGAGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAATGAGCTATTAGTCATGACGCTTTT");
        let mut aligner = SingleStrandAligner::default();
        aligner.scoring.jump_score = -1000;
        let alignment = aligner.global(&x, &y);
        let mut operations: Vec<AlignmentOperation> = Vec::new();
        operations.extend(&[Del; 19]);
        operations.extend(&[Match; 54]);
        operations.extend(&[Del; 27]);
        assert_alignment(
            &alignment,
            0,
            54,
            0,
            100,
            54 - (5 + (1 * 19)) - (5 + (1 * 27)),
            &operations,
        );
    }

    #[rstest]
    fn test_semiglobal_leading_and_trailing_deletions() {
        let x = s("-------------------GGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTG---------------------------");
        let y = s("AGGGCTATAGACTGCTAGAGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAATGAGCTATTAGTCATGACGCTTTT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.semiglobal(&x, &y);
        assert_alignment(&alignment, 0, 54, 19, 73, 54, &[Match; 54]);
    }

    #[rstest]
    fn test_local_identical() {
        let x = s("ACGTAACC");
        let y = s("ACGTAACC");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8, &[Match; 8]);
    }

    #[rstest]
    fn test_local_identical_query_in_target() {
        let x = s("  CCGG  ");
        let y = s("AACCGGTT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 0, 4, 2, 6, 4, &[Match; 4]);
    }

    #[rstest]
    fn test_local_identical_target_in_query() {
        let x = s("AACCGGTT");
        let y = s("  CCGG  ");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 2, 6, 0, 4, 4, &[Match; 4]);
    }

    #[rstest]
    fn test_local_y_subsequence_with_a_leading_mismatch() {
        // NB: first mismatch is not aligned
        let x = s("AGCGTCGTATACGTCGTA       ");
        let y = s("cGCGTCGTATACGTCGTAAAGATAT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 1, 18, 1, 18, 17, &[Match; 17]);
    }

    #[rstest]
    fn test_local_y_subsequence_with_a_trailing_mismatch() {
        let x = s("       CGCGTCGTATACGTCGTT");
        let y = s("AAGATATCGCGTCGTATACGTCGTa");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 0, 17, 7, 24, 17, &[Match; 17]);
    }

    #[rstest]
    fn test_local_y_subsequence_with_a_gap_in_x() {
        let x = s("  CCGCG-CGCGC  ");
        let y = s("AACCGCGACGCGCTT");
        let mut aligner = SingleStrandAligner::default();
        aligner.scoring.gap_open = -3;
        let alignment = aligner.local(&x, &y);
        assert_alignment(
            &alignment,
            0,
            10,
            2,
            13,
            10 - (3 + 1),
            &[
                Match, Match, Match, Match, Match, Del, Match, Match, Match, Match, Match,
            ],
        );
    }

    #[rstest]
    fn test_local_y_subsequence_with_a_gap_in_y() {
        let x = s("AACCGCGACGCGCTT");
        let y = s("  CCGCG-CGCGC  ");
        let mut aligner = SingleStrandAligner::default();
        aligner.scoring.gap_open = -3;
        let alignment = aligner.local(&x, &y);
        assert_alignment(
            &alignment,
            2,
            13,
            0,
            10,
            10 - (3 + 1),
            &[
                Match, Match, Match, Match, Match, Ins, Match, Match, Match, Match, Match,
            ],
        );
    }

    #[rstest]
    fn test_local_prefer_match_over_indel() {
        let x = s("  CGCGCGCG   ");
        //                    ||||
        let y = s("AACGCGACGCGTT");
        let mut aligner: SingleStrandAligner<MatchParams> = SingleStrandAligner::default();
        aligner.scoring.gap_open = -3;
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 0, 4, 2, 6, 4, &[Match; 4]);
    }

    #[rstest]
    fn test_local_zero_length_alignment() {
        let x = s("TTTTT");
        let y = s("AAAAA");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 1, 1, 5, 5, 0, &[]);
    }

    #[rstest]
    fn test_global_jump_with_leading_and_trailing_matches() {
        // The first 13bp of x and y align, then last 13bp of x and y align.  The "GATCGATC"
        // subsequence in x is aligned twice:
        // y:  TTTTTGATCGATC     ==>     GATCGATCTTTTT
        // x: [TTTTTGATCGATCTTTTT] [TTTTTGATCGATCTTTTT]
        // y:
        let x = s("TTTTT________GATCGATCTTTTT");
        let y = s("TTTTTGATCGATCGATCGATCTTTTT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        let mut operations: Vec<AlignmentOperation> = Vec::new();
        operations.extend(&[Match; 13]);
        operations.push(Xskip(13));
        operations.extend(&[Match; 13]);
        assert_alignment(&alignment, 0, 18, 0, 26, 26 - 10, &operations);
    }

    #[rstest]
    fn test_semiglobal_jump_with_leading_and_trailing_matches() {
        // The first 13bp of x and y align, then last 13bp of x and y align.  The "GATCGATC"
        // subsequence in x is aligned twice:
        // x: [TTTTTGATCGATCTTTTT] [TTTTTGATCGATCTTTTT]
        // y:  TTTTTGATCGATC     ==>     GATCGATCTTTTT
        let x = s("TTTTT________GATCGATCTTTTT");
        let y = s("TTTTTGATCGATCGATCGATCTTTTT");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.semiglobal(&x, &y);
        let mut operations: Vec<AlignmentOperation> = Vec::new();
        operations.extend(&[Match; 13]);
        operations.push(Xskip(13));
        operations.extend(&[Match; 13]);
        assert_alignment(&alignment, 0, 18, 0, 26, 26 - 10, &operations);
    }

    #[rstest]
    fn test_global_jump_back_to_start_of_x() {
        // The first 13bp of x and y align, then last 13bp of x and y align.  The "GATCGATC"
        // subsequence in x is aligned twice:
        // x: [GATCGATC] [GATCGATC]
        // y:  GATCGATC==>GATCGATC
        let x = s("GATCGATC________");
        let y = s("GATCGATCGATCGATC");
        let mut aligner = SingleStrandAligner::default();
        let alignment = aligner.global(&x, &y);
        let mut operations: Vec<AlignmentOperation> = Vec::new();
        operations.extend(&[Match; 8]);
        operations.push(Xskip(8));
        operations.extend(&[Match; 8]);
        assert_alignment(&alignment, 0, 8, 0, 16, 16 - 10, &operations);
    }
}
