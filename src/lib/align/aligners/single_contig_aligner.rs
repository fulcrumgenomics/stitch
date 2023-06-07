// Original Code was copied from:
//     https://github.com/rust-bio/rust-bio/blob/master/src/alignment/pairwise/mod.rs
// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp::max;
use std::i32;
use std::iter::repeat;

use crate::align::aligners::constants::AlignmentMode;
use crate::align::scoring::Scoring;
use crate::align::traceback::TB_XJUMP;
use bio::alignment::pairwise::MatchFunc;
use bio::alignment::pairwise::MatchParams;
use bio::utils::TextSlice;

use crate::align::aligners::constants::AlignmentOperation;
use crate::align::aligners::constants::DEFAULT_ALIGNER_CAPACITY;
use crate::align::alignment::Alignment;

use super::constants::MIN_SCORE;
use super::JumpInfo;
use crate::align::traceback::traceback;
use crate::align::traceback::Cell;
use crate::align::traceback::Traceback;
use crate::align::traceback::TracebackCell;
use crate::align::traceback::TB_DEL;
use crate::align::traceback::TB_INS;
use crate::align::traceback::TB_MATCH;
use crate::align::traceback::TB_START;
use crate::align::traceback::TB_SUBST;
use crate::align::traceback::TB_XCLIP_PREFIX;
use crate::align::traceback::TB_XCLIP_SUFFIX;
use crate::align::traceback::TB_YCLIP_PREFIX;
use crate::align::traceback::TB_YCLIP_SUFFIX;

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
pub struct SingleContigAligner<F: MatchFunc> {
    pub I: [Vec<i32>; 2],
    pub D: [Vec<i32>; 2],
    pub S: [Vec<i32>; 2],
    pub Lx: Vec<usize>,
    pub Ly: Vec<usize>,
    pub Sn: Vec<i32>,
    pub traceback: Traceback,
    pub scoring: Scoring<F>,
    pub contig_idx: u32,
    pub circular: bool,
}

impl Default for SingleContigAligner<MatchParams> {
    fn default() -> Self {
        let match_fn = MatchParams::new(1, -1);
        SingleContigAligner::new(-5, -1, -10, match_fn)
    }
}

impl<F: MatchFunc> SingleContigAligner<F> {
    pub fn set_contig_idx(&mut self, contig_idx: usize) {
        self.contig_idx = contig_idx as u32;
    }

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
                let mut tb = Cell::default();
                tb.set_all(TB_START, 0);
                tb.set_s_all(TB_START, 0, self.contig_idx, 0);
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
                let mut tb = Cell::default();
                tb.set_all(TB_START, 0);
                tb.set_s_all(TB_START, 0, self.contig_idx, 0);
                if i == 1 {
                    self.I[k][i] = self.scoring.gap_open + self.scoring.gap_extend;
                    tb.set_i(TB_START, 1);
                } else {
                    // Insert all i characters
                    // Could either be a single long-insertion, or x-clipping then an insertion start
                    let i_score = self.scoring.gap_open + self.scoring.gap_extend * (i as i32);
                    let c_score =
                        self.scoring.xclip_prefix + self.scoring.gap_open + self.scoring.gap_extend; // Clip then insert
                    if i_score > c_score {
                        self.I[k][i] = i_score;
                        tb.set_i(TB_INS, i as u32);
                    } else {
                        self.I[k][i] = c_score;
                        tb.set_i(TB_XCLIP_PREFIX, 0);
                    }
                }

                if i == m {
                    // tb.set_s(TB_XCLIP_SUFFIX, i as u32);
                    tb.set_s(TB_XCLIP_SUFFIX, 0);
                } else {
                    self.S[k][i] = MIN_SCORE;
                }

                if self.I[k][i] > self.S[k][i] {
                    self.S[k][i] = self.I[k][i];
                    tb.set_s(TB_INS, i as u32);
                }

                if self.scoring.xclip_prefix > self.S[k][i] {
                    self.S[k][i] = self.scoring.xclip_prefix;
                    // tb.set_s(TB_XCLIP_PREFIX, i as u32);
                    tb.set_s(TB_XCLIP_PREFIX, 0);
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
        let mut tb = Cell::default();
        tb.set_s_all(TB_START, 0, self.contig_idx, 0);
        self.I[curr][0] = MIN_SCORE;

        // deletion
        if j == 1 {
            // deletion start
            self.D[curr][0] = self.scoring.gap_open + self.scoring.gap_extend;
            tb.set_d(TB_START, 1);
        } else {
            // Delete all j characters
            // Could either be a single long-deletion, or y-clipping then an insertion start
            let d_score = self.scoring.gap_open + self.scoring.gap_extend * (j as i32);
            let c_score =
                self.scoring.yclip_prefix + self.scoring.gap_open + self.scoring.gap_extend;
            if d_score > c_score {
                self.D[curr][0] = d_score;
                tb.set_d(TB_DEL, j as u32);
            } else {
                self.D[curr][0] = c_score;
                // tb.set_d(TB_YCLIP_PREFIX, j as u32);
                tb.set_d(TB_YCLIP_PREFIX, 0);
            }
        }
        if self.D[curr][0] > self.scoring.yclip_prefix {
            self.S[curr][0] = self.D[curr][0];
            tb.set_s(TB_DEL, j as u32);
        } else {
            self.S[curr][0] = self.scoring.yclip_prefix;
            // tb.set_s(TB_YCLIP_PREFIX, j as u32);
            tb.set_s(TB_YCLIP_PREFIX, 0);
        }

        // Track the score if we do suffix clip (y) from here
        if j == n && self.Sn[0] > self.S[curr][0] {
            self.S[curr][0] = self.Sn[0];
            // tb.set_s(TB_YCLIP_SUFFIX, (n + m) as u32);
            tb.set_s(TB_YCLIP_SUFFIX, 0);
        } else if self.S[curr][0] + self.scoring.yclip_suffix > self.Sn[0] {
            self.Sn[0] = self.S[curr][0] + self.scoring.yclip_suffix;
            self.Ly[0] = n - j;
        }

        self.traceback.set(0, j, tb);

        // Handle i > 0 cases
        for i in 1..=m {
            self.S[curr][i] = MIN_SCORE;
        }
    }

    /// Gets the jump score for a given cell in the matrix.
    fn get_jump_score_and_len(
        &self,
        m: usize,
        i: usize,
        j: usize,
        prev: usize,
        addend: i32,
        jump_info: JumpInfo,
    ) -> JumpInfo {
        // add the specific addend!
        let jump_info = {
            let mut info = jump_info;
            info.score += addend;
            info
        };

        // DO NOT consider a circular no-cost jump from the end (previous) to the start (current)
        if !self.circular || i != 1 {
            return jump_info;
        }

        // Do not jump from an Xclip
        let jump_from_end_tb = self.traceback.get(m, j - 1).get_s().tb;
        if jump_from_end_tb == TB_XCLIP_SUFFIX {
            return jump_info;
        }

        // Get the score of jumping from the end of the previous column to the start of the current
        // column
        let jump_from_end_score = self.S[prev][m] + addend;
        if jump_info.score > jump_from_end_score {
            return jump_info;
        }

        // If equal, tie-break on the longest alignment length
        let jump_from_end_s = self.traceback.get(m, j - 1).get_s();
        let jump_from_end_len = jump_from_end_s.len + 1;
        if jump_from_end_score == jump_info.score && jump_from_end_len <= jump_info.len {
            return jump_info;
        }

        // return the zero-cost jump from the end
        JumpInfo {
            score: jump_from_end_score,
            len: jump_from_end_len,
            idx: self.contig_idx,
            from: m as u32,
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
        jump_info: JumpInfo,
    ) {
        let q = y[j - 1];
        let xclip_score = self.scoring.xclip_prefix
            + max(
                self.scoring.yclip_prefix,
                self.scoring.gap_open + self.scoring.gap_extend * (j as i32),
            );

        for i in 1..=m {
            let p: u8 = x[i - 1];
            let mut tb = Cell::default();

            // Insertion
            // It does not make sense to _start_ an insertion right after a jump, since you might
            // as well just jumped over the insertion!
            let i_score = self.I[curr][i - 1] + self.scoring.gap_extend;
            let s_score: i32 =
                self.S[curr][i - 1] + self.scoring.gap_open + self.scoring.gap_extend;
            let best_i_score = max(i_score, s_score);
            if i_score == best_i_score {
                tb.set_i(TB_INS, self.traceback.get(i - 1, j).get_i_len() + 1);
            } else {
                let s_value = self.traceback.get(i - 1, j).get_s();
                tb.set_i(s_value.tb, s_value.len + 1);
            }

            // Deletion
            let d_score = self.D[prev][i] + self.scoring.gap_extend;
            let s_score = self.S[prev][i] + self.scoring.gap_open + self.scoring.gap_extend;
            let best_d_score = max(d_score, s_score);
            if d_score == best_d_score {
                let prev_len = self.traceback.get(i, j - 1).get_d_len();
                tb.set_d(TB_DEL, prev_len + 1);
            } else {
                let s_value = self.traceback.get(i, j - 1).get_s();
                tb.set_d(s_value.tb, s_value.len + 1);
            }

            // Set the optimal score for all moves
            // Preferences if two or more moves have
            // 1. diagonal over all other moves except jump, and over jump when the alignment length
            // for the diagonal is greater than the jump's alignment length.
            // 2. X-suffix clip (for implementation convenience)
            // 3. deletion
            // 4. insertion
            // 5. jump *(exception see rule 1)
            // 6. X-prefix clip
            // 7. Y-prefix clip
            tb.set_s(TB_XCLIP_SUFFIX, self.traceback.get(i, j).get_s_len());
            let mut best_s_score = self.S[curr][i];
            // Score for aligning just [x-1] with y[j-1] alone
            let addend = self.scoring.match_fn.score(p, q);
            // Align the x[i-1] with y[j-1] through a diagonal move.
            let diag_score = self.S[prev][i - 1] + addend;
            let diag_len = self.traceback.get(i - 1, j - 1).get_s_len() + 1;
            if diag_score >= best_s_score {
                best_s_score = diag_score;
                let s_tb = if p == q { TB_MATCH } else { TB_SUBST };
                tb.set_s_all(s_tb, diag_len, self.contig_idx, (i - 1) as u32);
            }
            // Deletion
            if best_d_score > best_s_score {
                best_s_score = best_d_score;
                tb.set_s_all(TB_DEL, tb.get_d_len(), self.contig_idx, i as u32);
            }
            // Insertion
            if best_i_score > best_s_score {
                best_s_score = best_i_score;
                tb.set_s_all(TB_INS, tb.get_i_len(), self.contig_idx, (i - 1) as u32);
            }
            // Align the x[i-1] with y[j-1] through a jump move.
            let x_jump_info = self.get_jump_score_and_len(m, i, j, prev, addend, jump_info);
            let do_jump = x_jump_info.score > best_s_score
                || (x_jump_info.score == best_s_score
                    && best_s_score == diag_score
                    && x_jump_info.len > diag_len);
            if do_jump {
                best_s_score = x_jump_info.score;
                let s_tb = if p == q { TB_MATCH } else { TB_SUBST };
                tb.set_s_all(s_tb, x_jump_info.len, x_jump_info.idx, x_jump_info.from);
            }
            // X-prefix clip
            if xclip_score > best_s_score {
                best_s_score = xclip_score;
                let prev_len = self.traceback.get(0, j).get_s_len();
                // tb.set_s_all(TB_XCLIP_PREFIX, prev_len + i as u32, 0, false);
                tb.set_s_all(TB_XCLIP_PREFIX, prev_len, self.contig_idx, 0);
            }
            // Y-prefix clip
            let yclip_score = self.scoring.yclip_prefix
                + self.scoring.gap_open
                + self.scoring.gap_extend * (i as i32);
            if yclip_score > best_s_score {
                let prev_len = self.traceback.get(i, 0).get_s_len();
                best_s_score = yclip_score;
                // tb.set_s_all(TB_YCLIP_PREFIX, prev_len + j as u32, i as u32, false);
                tb.set_s_all(TB_YCLIP_PREFIX, prev_len, self.contig_idx, i as u32);
            }

            // Set the values in the matrices
            self.S[curr][i] = best_s_score;
            self.I[curr][i] = best_i_score;
            self.D[curr][i] = best_d_score;

            // Track the score if we do suffix clip (x) from here
            let do_x_suffix_clip =
                match (self.S[curr][i] + self.scoring.xclip_suffix).cmp(&self.S[curr][m]) {
                    std::cmp::Ordering::Less => false,
                    std::cmp::Ordering::Greater => true,
                    std::cmp::Ordering::Equal => {
                        // let left_len = tb.get_s_len() + (m - i) as u32;
                        let left_len = tb.get_s_len();
                        let right_len = self.traceback.get(m, j).get_s_len();
                        left_len > right_len
                    }
                };
            if do_x_suffix_clip {
                self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
                let prev_s: crate::align::traceback::SValue = tb.get_s();
                self.traceback.get_mut(m, j).set_s_all(
                    TB_XCLIP_SUFFIX,
                    // prev_len + (m - i) as u32,
                    prev_s.len,
                    prev_s.idx,
                    i as u32,
                );
                self.Lx[j] = m - i;
            }

            // Track the score if we do suffix clip (y) from here
            let do_y_suffix_clip =
                match (self.S[curr][i] + self.scoring.yclip_suffix).cmp(&self.Sn[i]) {
                    std::cmp::Ordering::Less => false,
                    std::cmp::Ordering::Greater => true,
                    std::cmp::Ordering::Equal => {
                        // let left_len = tb.get_s_len() + (n - j) as u32;
                        let left_len = tb.get_s_len();
                        let right_len = self.traceback.get(i, n).get_s_len();
                        left_len > right_len
                    }
                };

            if do_y_suffix_clip {
                self.Sn[i] = self.S[curr][i] + self.scoring.yclip_suffix;
                self.Ly[i] = n - j;
            }

            self.traceback.set(i, j, tb);
        }
    }

    pub fn fill_last_column_and_end_clipping(&mut self, m: usize, n: usize) {
        // Handle jumping over the remaining i bases in x and suffix clipping, in the j=n case
        for i in 0..=m {
            let j: usize = n; // end of y
            let curr: usize = j % 2;

            // jump over the remaining i bases in x
            if self.S[curr][i] + self.scoring.jump_score_same_contig_and_strand > self.S[curr][m] {
                self.S[curr][m] = self.S[curr][i] + self.scoring.jump_score_same_contig_and_strand;
                let prev_s = self.traceback.get(i, j).get_s();
                self.traceback
                    .get_mut(m, j)
                    .set_s_all(TB_XJUMP, prev_s.len, prev_s.idx, i as u32);
            }

            // y-clip
            let do_y_suffix_clip = match (self.Sn[i]).cmp(&self.S[curr][i]) {
                std::cmp::Ordering::Less => false,
                std::cmp::Ordering::Greater => true,
                std::cmp::Ordering::Equal => {
                    let left_len = self.traceback.get(i, n).get_s_len();
                    // let right_len: u32 = self.traceback.get(i, j).get_s_len() + (n - j) as u32;
                    let right_len: u32 = self.traceback.get(i, j).get_s_len();
                    left_len > right_len
                }
            };
            if do_y_suffix_clip {
                self.S[curr][i] = self.Sn[i];
                // no need to set Ly[i] since it's already set in fill_last_column
                let s_value = self.traceback.get(i, j - self.Ly[i]).get_s();
                let tb = self.traceback.get_mut(i, j);
                tb.set_s_all(
                    TB_YCLIP_SUFFIX,
                    // s_value.len + self.Ly[i] as u32,
                    s_value.len,
                    s_value.idx,
                    i as u32,
                );
            }

            // x-clip
            let do_x_suffix_clip =
                match (self.S[curr][i] + self.scoring.xclip_suffix).cmp(&self.S[curr][m]) {
                    std::cmp::Ordering::Less => false,
                    std::cmp::Ordering::Greater => true,
                    std::cmp::Ordering::Equal => {
                        // let left_len = self.traceback.get(i, j).get_s_len() + (m - i) as u32;
                        let left_len = self.traceback.get(i, j).get_s_len();
                        let right_len = self.traceback.get(m, j).get_s_len();
                        left_len > right_len
                    }
                };
            if do_x_suffix_clip {
                self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
                self.Lx[j] = m - i;
                let prev_s = self.traceback.get(i, j).get_s();
                self.traceback.get_mut(m, j).set_s_all(
                    TB_XCLIP_SUFFIX,
                    // prev_len + (m - i) as u32,
                    prev_s.len,
                    prev_s.idx,
                    i as u32,
                );
            }
        }

        // Since there could be a change in the last column of S,
        // recompute the last column of I as this could also change
        for i in 1..=m {
            let j = n;
            let curr = j % 2;
            let i_score = self.S[curr][i - 1] + self.scoring.gap_open + self.scoring.gap_extend;
            if i_score > self.I[curr][i] {
                self.I[curr][i] = i_score;
                let s_value = self.traceback.get(i - 1, j).get_s();
                self.traceback
                    .get_mut(i, j)
                    .set_i(s_value.tb, s_value.len + 1);
            }

            if i_score > self.S[curr][i] {
                self.S[curr][i] = i_score;
                let prev_len = self.traceback.get_mut(i, j).get_i_len();
                self.traceback.get_mut(i, j).set_s_all(
                    TB_INS,
                    prev_len,
                    self.contig_idx,
                    (i - 1) as u32,
                );
                if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][m] {
                    self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
                    self.Lx[j] = m - i;
                    self.traceback.get_mut(m, j).set_s_all(
                        TB_XCLIP_SUFFIX,
                        // prev_len + ((m - i) as u32),
                        prev_len,
                        self.contig_idx,
                        i as u32,
                    );
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
        SingleContigAligner::with_capacity(
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

        SingleContigAligner {
            I: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            D: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            S: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            Lx: Vec::with_capacity(n + 1),
            Ly: Vec::with_capacity(m + 1),
            Sn: Vec::with_capacity(m + 1),
            traceback: Traceback::with_capacity(m, n),
            scoring: Scoring::with_jump_score(gap_open, gap_extend, jump_score, match_fn),
            contig_idx: 0,
            circular: false,
        }
    }

    /// Create new aligner instance with given the scoring struct
    ///
    /// # Arguments
    ///
    /// * `scoring` - the scoring struct (see bio::alignment::pairwise::Scoring)
    #[allow(dead_code)]
    pub fn with_scoring(scoring: Scoring<F>) -> Self {
        SingleContigAligner::with_capacity_and_scoring(
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

        SingleContigAligner {
            I: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            D: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            S: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            Lx: Vec::with_capacity(n + 1),
            Ly: Vec::with_capacity(m + 1),
            Sn: Vec::with_capacity(m + 1),
            traceback: Traceback::with_capacity(m, n),
            scoring,
            contig_idx: 0,
            circular: false,
        }
    }

    /// Sets the value for treating x as circular, allowing for a zero-cost jump to the start of x.
    pub fn set_circular(&mut self, circular: bool) {
        self.circular = circular;
    }

    /// Gets the best jump score and x-index for the jump
    pub fn get_jump_info(&self, m: usize, j: usize, jump_score: i32) -> JumpInfo {
        let cur = j % 2;

        let mut best_jump_score = self.S[cur][0] + jump_score;
        let mut best_jump_from = 0;
        for k in 1..=m {
            if best_jump_score < self.S[cur][k] + jump_score {
                best_jump_score = self.S[cur][k] + jump_score;
                best_jump_from = k;
            }
        }

        let best_jump_len = self.traceback.get(best_jump_from, j).get_s_len() + 1;

        JumpInfo {
            score: best_jump_score,
            from: best_jump_from as u32,
            idx: self.contig_idx,
            len: best_jump_len,
        }
    }

    /// The core function to compute the alignment
    ///
    /// # Arguments
    ///
    /// * `x` - Textslice
    /// * `y` - Textslice
    pub fn custom(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        let (m, n) = (x.len(), y.len());

        self.init_matrices(m, n);

        for j in 1..=n {
            let curr = j % 2;
            let prev = 1 - curr;

            // Initialize the column
            self.init_column(j, curr, m, n);

            // Get the best jump score and x-index for the jump
            let jump_info =
                self.get_jump_info(m, j - 1, self.scoring.jump_score_same_contig_and_strand);

            // Fill the column
            self.fill_column(x, y, m, n, j, prev, curr, jump_info);
        }

        self.fill_last_column_and_end_clipping(m, n);

        let aligners = vec![&*self];
        traceback(&aligners, n)
    }

    /// Calculate global alignment of x against y.
    #[allow(dead_code)]
    pub fn global(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
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
    #[allow(dead_code)]
    pub fn querylocal(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
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
        alignment.mode = AlignmentMode::QueryLocal;

        // Filter out Yclip from alignment.operations
        {
            use self::AlignmentOperation::Yclip;
            alignment.operations.retain(|x| !matches!(*x, Yclip(_)));
        }

        // Set the clip penalties to the original values
        self.scoring.xclip_prefix = clip_penalties[0];
        self.scoring.xclip_suffix = clip_penalties[1];
        self.scoring.yclip_prefix = clip_penalties[2];
        self.scoring.yclip_suffix = clip_penalties[3];

        alignment
    }

    /// Calculate semiglobal alignment of x against y (x is local, y is global).
    #[allow(dead_code)]
    pub fn targetlocal(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
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
        self.scoring.yclip_prefix = MIN_SCORE;
        self.scoring.yclip_suffix = MIN_SCORE;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::TargetLocal;

        // Filter out Xclip from alignment.operations
        {
            use self::AlignmentOperation::Xclip;
            alignment.operations.retain(|x| !matches!(*x, Xclip(_)));
        }

        // Set the clip penalties to the original values
        self.scoring.xclip_prefix = clip_penalties[0];
        self.scoring.xclip_suffix = clip_penalties[1];
        self.scoring.yclip_prefix = clip_penalties[2];
        self.scoring.yclip_suffix = clip_penalties[3];

        alignment
    }

    /// Calculate local alignment of x against y.
    #[allow(dead_code)]
    pub fn local(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
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
            // Filter out Xclip from alignment.operations
            {
                use self::AlignmentOperation::{Xclip, Yclip};
                alignment
                    .operations
                    .retain(|x| !matches!(*x, Xclip(_) | Yclip(_)));
            }
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
    use bio::alignment::pairwise::MatchParams;
    use itertools::Itertools;
    use rstest::rstest;

    use crate::align::alignment::Alignment;

    use super::SingleContigAligner;

    /// Upper-cases and remove display-related characters from a string.
    fn s(bases: &str) -> Vec<u8> {
        bases
            .chars()
            .filter(|base| *base != '-' && *base != ' ' && *base != '_')
            .map(|base| base.to_ascii_uppercase() as u8)
            .collect_vec()
    }

    fn assert_alignment(
        alignment: &Alignment,
        xstart: usize,
        xend: usize,
        ystart: usize,
        yend: usize,
        score: i32,
        cigar: &str,
        length: usize,
    ) {
        assert_eq!(alignment.xstart, xstart, "xstart {alignment}");
        assert_eq!(alignment.xend, xend, "xend {alignment}");
        assert_eq!(alignment.ystart, ystart, "ystart {alignment}");
        assert_eq!(alignment.yend, yend, "yend {alignment}");
        assert_eq!(alignment.score, score, "score {alignment}");
        assert_eq!(alignment.start_contig_idx, 0, "strand {alignment}");
        assert_eq!(alignment.cigar(), cigar, "cigar {alignment}");
        assert_eq!(alignment.length, length, "length {alignment}");
    }

    /// Identical sequences, all matches
    #[rstest]
    fn test_identical() {
        let x = s("ACGTAACC");
        let y = s("ACGTAACC");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8, "8=", 8);
    }

    /// Identical sequences, except one mismatch
    #[rstest]
    fn test_single_mismatch() {
        let x = s("AACCGGTT");
        let y = s("AACCGtTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 8, 0, 8, 7 - 1, "5=1X2=", 8);
    }

    /// Identical sequences, except one samll deletion
    #[rstest]
    fn test_small_deletion() {
        let x = s("AACC-GTT");
        let y = s("AACCGGTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 7, 0, 8, 7 - (5 + 1), "4=1D3=", 8);
    }

    /// Identical sequences, except one samll insertion
    #[rstest]
    fn test_small_insertion() {
        let x = s("AACCGGTT");
        let y = s("AACC-GTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 8, 0, 7, 7 - (5 + 1), "4=1I3=", 8);
    }

    /// Two sequences with compensating insertions and deletions
    #[rstest]
    fn test_compensating_insertion_and_deletion() {
        let x = s("AAACGCGCGCGCG-TT");
        let y = s("-AACGCGCGCGCGTTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            15,
            0,
            15,
            14 - (5 + 1) - (5 + 1),
            "1I12=1D2=",
            16,
        );
    }

    #[rstest]
    fn test_leading_insertion() {
        let x = s("ATTTTTTTTTTT");
        let y = s("-TTTTTTTTTTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 12, 0, 11, 11 - (5 + 1), "1I11=", 12);
    }

    #[rstest]
    fn test_trailing_insertion() {
        let x = s("TTTTTTTTTTTA");
        let y = s("TTTTTTTTTTT-");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 12, 0, 11, 11 - (5 + 1), "11=1I", 12);
    }

    #[rstest]
    fn test_leading_deletion() {
        let x = s("-TTTTTTTTTTT");
        let y = s("ATTTTTTTTTTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 11, 0, 12, 11 - (5 + 1), "1D11=", 12);
    }

    #[rstest]
    fn test_trailing_deletion() {
        let x = s("TTTTTTTTTTT-");
        let y = s("TTTTTTTTTTTA");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 11, 0, 12, 11 - (5 + 1), "11=1D", 12);
    }

    /// Align two sequences preferring a 2bp insertion and mismatch vs two small insertions"
    #[rstest]
    fn test_prefer_2bp_insertion_and_mismatch_vs_two_small_insertions() {
        let x = s("ATTTTTTTTTTTA");
        let y = s("--TTTTTTTTTTt");
        let match_fn = MatchParams::new(1, -1);
        let mut aligner = SingleContigAligner::new(-3, -1, -10, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            13,
            0,
            11,
            10 - (3 + 1) - 1 - 1,
            "2I10=1X",
            13,
        );
    }

    /// Align two sequences preferring a 2bp insertion and mismatch vs two small insertions"
    #[rstest]
    fn test_prefer_two_small_insertions_vs_2bp_insertion_and_mismatch() {
        let x = s("ATTTTTTTTTTTA");
        let y = s("-TTTTTTTTTTT-");
        let match_fn = MatchParams::new(1, -3);
        let mut aligner = SingleContigAligner::new(-3, -1, -10, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            13,
            0,
            11,
            11 - (3 + 1) - (3 + 1),
            "1I11=1I",
            13,
        );
    }

    #[rstest]
    fn test_left_justify_insertion_in_homopolymer() {
        let x = s("GTTTTTTTTTTA");
        let y = s("G-TTTTTTTTTA");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 12, 0, 11, 11 - (5 + 1), "1=1I10=", 12);
    }

    #[rstest]
    fn test_left_justify_insertion_in_triplet() {
        let x = s("GACGACGACGACGA");
        let y = s("---GACGACGACGA");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 14, 0, 11, 11 - (5 + 1) - 1 - 1, "3I11=", 14);
    }

    #[rstest]
    fn test_left_justify_insertion_in_triplet_with_leading_matches() {
        let x = s("TTTGACGACGACGACGA");
        let y = s("TTT---GACGACGACGA");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            17,
            0,
            14,
            14 - (5 + 1) - 1 - 1,
            "3=3I11=",
            17,
        );
    }

    /// NB: if the jump score is set to -11, then the alignment would jump back in x (3bp), to give
    /// 17 matches (+17) and one jump (-11) for a score of 6, which is equal to 14 matches (+14)
    /// and a 3bp deletion (-5 - 1 - 1 -1 = -8) for a score of 6.
    #[rstest]
    fn test_jump_over_deletion_in_triplet() {
        let x = s("TTTGACGACGA___CGA");
        let y = s("TTTGACGACGACGACGA");
        let mut aligner = SingleContigAligner::default();
        aligner.scoring = aligner.scoring.set_jump_score(-11);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            14,
            0,
            17,
            14 - (5 + 1 + 1 + 1),
            "3=3D11=",
            17,
        );
    }

    /// NB: if the jump score is set to -11, then the alignment would jump back in x (3bp), to give
    /// 17 matches (+17) and one jump (-11) for a score of 6, which is equal to 14 matches (+14)
    /// and a 3bp deletion (-5 - 1 - 1 -1 = -8) for a score of 6.  We prefer the deletion in this case.
    #[rstest]
    fn test_deletion_over_jump() {
        let x = s("TTT---GACGACGACGA");
        let y = s("TTTGACGACGACGACGA");
        let mut aligner = SingleContigAligner::default();
        aligner.scoring = aligner.scoring.set_jump_score(-11);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            14,
            0,
            17,
            14 - (5 + 1 + 1 + 1),
            "3=3D11=",
            17,
        );
    }

    /// NB: if the jump score is set to -10, then the alignment would jump back in x (3bp), to give
    /// 17 matches (+17) and one jump (-10) for a score of 7, which is greater than 14 matches (+14)
    /// and a 3bp deletion (-5 - 1 - 1 -1 = -8) for a score of 6.
    #[rstest]
    fn test_jump_over_deletion() {
        let x = s("TTT___GACGACGACGA");
        let y = s("TTTGACGACGACGACGA");
        let mut aligner = SingleContigAligner::default();
        aligner.scoring = aligner.scoring.set_jump_score(-10);
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 14, 0, 17, 17 - 10, "6=3j11=", 17);
    }

    /// Prefer a mismatch over an insertion and deletion when the mismatch score is less than than
    /// two times the gap open + gap extend scores. NB: could be either "1I2=1D3=" or "2=1X3=".
    #[rstest]
    fn test_prefer_mismatch_over_insertion_and_deletion() {
        let x = s("AAACCC");
        let y = s("AAcCCC");
        let match_fn = MatchParams::new(1, -3);
        let mut aligner = SingleContigAligner::new(-1, -1, -10, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 6, 0, 6, 5 - 3, "2=1X3=", 6);
    }

    #[rstest]
    fn test_prefer_mismatch_over_insertion_and_deletion_when_same_score() {
        let x = s("AAACCC");
        let y = s("AAcCCC");
        // NB: could be either "1I2=1D3=" or "2=1X3="
        let match_fn = MatchParams::new(1, -4);
        let mut aligner = SingleContigAligner::new(-1, -1, -10, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 6, 0, 6, 5 - 4, "2=1X3=", 6);
    }

    /// Prefer an insertion and deletion over a mismatch when the mismatch score is greater than
    /// two times the gap open + gap extend scores.
    #[rstest]
    fn test_prefer_insertion_and_deletion_over_mismatch() {
        let x = s("AAA-CCC");
        let y = s("AA-CCCC");
        // NB: could be either "1I2=1D3=" or "2=1X3="
        // NB: if we prefer a insertion over a deletion, then it would be 2M1D1I3M
        let match_fn = MatchParams::new(1, -5);
        let mut aligner = SingleContigAligner::new(-1, -1, -10, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 6, 0, 6, 5 - (1 + 1) - (1 + 1), "1I2=1D3=", 7);
    }

    #[rstest]
    fn test_prefer_one_insertion_with_large_gap_open() {
        let x = s("ATTTTTTTTTTTA");
        let y = s("--TTTTTTTTTTt");
        let match_fn = MatchParams::new(1, -5);
        let mut aligner = SingleContigAligner::new(-100, -1, -10000, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            13,
            0,
            11,
            10 - (100 + 1) - 1 - 5,
            "2I10=1X",
            13,
        );
    }

    #[rstest]
    fn test_prefer_two_small_insertions_with_large_gap_extend() {
        let x = s("ATTTTTTTTTTTA");
        let y = s("-TTTTTTTTTTT-");
        let match_fn = MatchParams::new(1, -5);
        let mut aligner = SingleContigAligner::new(-1, -100, -10000, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            13,
            0,
            11,
            11 - (100 + 1) - (100 + 1),
            "1I11=1I",
            13,
        );
    }

    #[rstest]
    fn test_querylocal_identical() {
        let x = s("ACGTAACC");
        let y = s("ACGTAACC");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.querylocal(&x, &y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8, "8=", 8);
    }

    #[rstest]
    fn test_querylocal_identical_subsequence() {
        let x = s("  CCGG  ");
        let y = s("AACCGGTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.querylocal(&x, &y);
        assert_alignment(&alignment, 0, 4, 2, 6, 4, "4=", 4);
    }

    #[rstest]
    fn test_querylocal_subsequence_with_mismatch() {
        let x = s("       CGCGTCGTATACGTCGTT");
        let y = s("AAGATATCGCGTCGTATACGTCGTa");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.querylocal(&x, &y);
        assert_alignment(&alignment, 0, 18, 7, 25, 17 - 1, "17=1X", 18);
    }

    #[rstest]
    fn test_querylocal_subsequence_with_deletion() {
        let x = s("  CGCG-CGCG  ");
        let y = s("AACGCGACGCGTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.querylocal(&x, &y);
        assert_alignment(&alignment, 0, 8, 2, 11, 8 - (5 + 1), "4=1D4=", 9);
    }

    #[rstest]
    fn test_querylocal_insertion_when_x_longer_than_y() {
        let x = s("AAAAGGGGTTTT");
        let y = s("AAAA----TTTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.querylocal(&x, &y);
        assert_alignment(
            &alignment,
            0,
            12,
            0,
            8,
            8 - (5 + 1) - 1 - 1 - 1,
            "4=4I4=",
            12,
        );
    }

    #[rstest]
    fn test_global_leading_and_trailing_deletions() {
        let x = s("-------------------GGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTG---------------------------");
        let y = s("AGGGCTATAGACTGCTAGAGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAATGAGCTATTAGTCATGACGCTTTT");
        let mut aligner = SingleContigAligner::default();
        aligner.scoring = aligner.scoring.set_jump_score(-1000);
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            54,
            0,
            100,
            54 - (5 + (1 * 19)) - (5 + (1 * 27)),
            "19D54=27D",
            100,
        );
    }

    #[rstest]
    fn test_querylocal_leading_and_trailing_deletions() {
        let x = s("-------------------GGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTG---------------------------");
        let y = s("AGGGCTATAGACTGCTAGAGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAATGAGCTATTAGTCATGACGCTTTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.querylocal(&x, &y);
        assert_alignment(&alignment, 0, 54, 19, 73, 54, "54=", 54);
    }

    #[rstest]
    fn test_local_identical() {
        let x = s("ACGTAACC");
        let y = s("ACGTAACC");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8, "8=", 8);
    }

    #[rstest]
    fn test_local_identical_query_in_target() {
        let x = s("  CCGG  ");
        let y = s("AACCGGTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 0, 4, 2, 6, 4, "4=", 4);
    }

    #[rstest]
    fn test_local_identical_target_in_query() {
        let x = s("AACCGGTT");
        let y = s("  CCGG  ");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 2, 6, 0, 4, 4, "4=", 4);
    }

    #[rstest]
    fn test_local_y_subsequence_with_a_leading_mismatch() {
        // NB: first mismatch is not aligned
        let x = s("AGCGTCGTATACGTCGTA       ");
        let y = s("cGCGTCGTATACGTCGTAAAGATAT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 1, 18, 1, 18, 17, "17=", 17);
    }

    #[rstest]
    fn test_local_y_subsequence_with_a_trailing_mismatch() {
        let x = s("       CGCGTCGTATACGTCGTT");
        let y = s("AAGATATCGCGTCGTATACGTCGTa");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 0, 17, 7, 24, 17, "17=", 17); // NB: suffix clipping counts for two
    }

    #[rstest]
    fn test_local_y_subsequence_with_a_gap_in_x() {
        let x = s("  CCGCG-CGCGC  ");
        let y = s("AACCGCGACGCGCTT");
        let mut aligner = SingleContigAligner::default();
        aligner.scoring.gap_open = -3;
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 0, 10, 2, 13, 10 - (3 + 1), "5=1D5=", 11);
    }

    #[rstest]
    fn test_local_y_subsequence_with_a_gap_in_y() {
        let x = s("AACCGCGACGCGCTT");
        let y = s("  CCGCG-CGCGC  ");
        let mut aligner = SingleContigAligner::default();
        aligner.scoring.gap_open = -3;
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 2, 13, 0, 10, 10 - (3 + 1), "5=1I5=", 11);
    }

    #[rstest]
    fn test_local_prefer_match_over_indel() {
        let x = s("       CGCGCGCG");
        //                         ||||
        let y = s("AACGCGACGCGTT  ");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        aligner.scoring.gap_open = -3;
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 0, 4, 7, 11, 4, "4=", 4);
    }

    #[rstest]
    fn test_local_zero_length_alignment() {
        let x = s("TTTTT");
        let y = s("AAAAA");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 0, 0, 0, 0, 0, "", 0);
    }

    #[rstest]
    fn test_global_jump_with_leading_and_trailing_matches() {
        // The first 13bp of x and y align, then last 13bp of x and y align.  The "GATCGATC"
        // subsequence in x is aligned twice:
        //    [0 ....... 13]    [5 ...      17]
        // x: [TTTTTGATCGAT]    [GATCGATCTTTTT]
        // y:  TTTTTGATCGAT  ==> GATCGATCTTTTT
        let x = s("TTTTTGATCGAT________CTTTTT");
        let y = s("TTTTTGATCGATCGATCGATCTTTTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 18, 0, 26, 26 - 10, "13=8j13=", 26);
    }

    #[rstest]
    fn test_querylocal_jump_with_leading_and_trailing_matches() {
        // The first 13bp of x and y align, then last 13bp of x and y align.  The "GATCGATC"
        // subsequence in x is aligned twice:
        //    [0 ........ 12]   [5 ...      17]
        // x: [TTTTTGATCGATC]   [GATCGATCTTTTT]
        // y:  TTTTTGATCGATC ==> GATCGATCTTTTT
        let x = s("TTTTT________GATCGATCTTTTT");
        let y = s("TTTTTGATCGATCGATCGATCTTTTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.querylocal(&x, &y);

        assert_alignment(&alignment, 0, 18, 0, 26, 26 - 10, "13=8j13=", 26);
    }

    #[rstest]
    fn test_global_jump_back_to_start_of_x() {
        // The first 13bp of x and y align, then last 13bp of x and y align.  The "GATCGATC"
        // subsequence in x is aligned twice:
        //    [0 .... 7] [0 .... 7]
        // x: [GATCGATC] [GATCGATC]
        // y:  GATCGATC==>GATCGATC
        let x = s("GATCGATC________");
        let y = s("GATCGATCGATCGATC");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 8, 0, 16, 16 - 10, "8=8j8=", 16);
    }

    #[rstest]
    fn test_global_triple_jump_back_to_start_of_x() {
        // The first 13bp of x and y align, then last 13bp of x and y align.  The "GATCGATC"
        // subsequence in x is aligned twice:
        //    [0 .... 7] [0 .... 7] [0 .... 7]
        // x: [GATCGATC] [GATCGATC] [GATCGATC]
        // y:  GATCGATC==>GATCGATC==>GATCGATC
        let x = s("GATCGATC________________");
        let y = s("GATCGATCGATCGATCGATCGATC");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 8, 0, 24, 24 - 10 - 10, "8=8j8=8j8=", 24);
    }

    #[rstest]
    fn test_global_sir_jump_a_lot() {
        //    [0 ...... 9] [20 .... 29] [10 .... 19] [30 .... 39]
        // x: [AAAAAAAAAA] [CCCCCCCCCC] [GGGGGGGGGG] [TTTTTTTTTT]
        // y:  AAAAAAAAAA==>CCCCCCCCCC==>GGGGGGGGGG==>TTTTTTTTTT
        //    [0 ...... 9] [10 .... 19] [20 .... 29] [30 .... 39]
        let x = s("AAAAAAAAAAGGGGGGGGGGCCCCCCCCCCTTTTTTTTTT");
        let y = s("AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(
            &alignment,
            0,
            40,
            0,
            40,
            40 - 10 - 10 - 10,
            "10=10J10=20j10=10J10=",
            40,
        );
    }

    #[rstest]
    fn test_querylocal_sir_jump_a_lot() {
        //    [0 ...... 9] [20 .... 29] [10 .... 19] [30 .... 39]
        // x: [AAAAAAAAAA] [CCCCCCCCCC] [GGGGGGGGGG] [TTTTTTTTTT]
        // y:  AAAAAAAAAA==>CCCCCCCCCC==>GGGGGGGGGG==>TTTTTTTTTT
        //    [0 ...... 9] [10 .... 19] [20 .... 29] [30 .... 39]
        let x = s("AAAAAAAAAAGGGGGGGGGGCCCCCCCCCCTTTTTTTTTT");
        let y = s("AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT");
        let mut aligner = SingleContigAligner::default();
        let alignment = aligner.querylocal(&x, &y);
        assert_alignment(
            &alignment,
            0,
            40,
            0,
            40,
            40 - 10 - 10 - 10,
            "10=10J10=20j10=10J10=",
            40,
        );
    }

    #[rstest]
    fn test_local_sir_jump_a_lot() {
        // NB: prefer longer alignments with jumps versus shorter alignments with fewer/no jumps
        // when the score is tied
        //    [0 ...... 9] [20 .... 29] [10 .... 19] [30 .... 39]
        // x: [AAAAAAAAAA] [CCCCCCCCCC] [GGGGGGGGGG] [TTTTTTTTTT]
        // y:  AAAAAAAAAA==>CCCCCCCCCC==>GGGGGGGGGG==>TTTTTTTTTT
        //    [0 ...... 9] [10 .... 19] [20 .... 29] [30 .... 39]
        let x = s("AAAAAAAAAAGGGGGGGGGGCCCCCCCCCCTTTTTTTTTT");
        let y = s("AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(
            &alignment,
            0,
            40,
            0,
            40,
            40 - 10 - 10 - 10,
            "10=10J10=20j10=10J10=",
            40,
        );
    }

    #[rstest]
    fn test_local_prefer_suffix_clip_to_last_jump() {
        // NB: only 9 Cs match at the start of x and end of y, so not enough matches to jump after
        // matching the 10 As at the end of x and start of y
        //    [9 ..... 18]
        // x: [AAAAAAAAAA
        // y:  AAAAAAAAAA
        //    [0 ...... 9]
        let x = s("CCCCCCCCCAAAAAAAAAA");
        let y = s("AAAAAAAAAACCCCCCCCC");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 9, 19, 0, 10, 10, "10=", 10);
    }

    #[rstest]
    fn test_local_prefer_last_jump_to_suffix_clip() {
        // NB: 10 Cs match at the start of x and end of y, so enough matches to jump after
        // matching the 10 As at the end of x and start of y.
        //    [9 ..... 19] [0 ...... 9]
        // x: [AAAAAAAAAA] [CCCCCCCCCC]
        // y:  AAAAAAAAAA==>CCCCCCCCCC
        //    [0 ...... 9] [10 .... 19]
        let x = s("CCCCCCCCCCAAAAAAAAAA");
        let y = s("AAAAAAAAAACCCCCCCCCC");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 10, 10, 0, 20, 20 - 10, "10=20j10=", 20);
    }

    #[rstest]
    fn test_local_prefer_prefix_clip_to_last_jump() {
        // NB: only 9 Cs match at the end of x and start of y, so not enough matches to jump after
        // matching the 10 As at the start of x and end of x.
        //    [0 ...... 9]
        // x: [AAAAAAAAAA
        // y:  AAAAAAAAAA
        //    [9 ...... 18]
        let x = s("AAAAAAAAAACCCCCCCCC");
        let y = s("CCCCCCCCCAAAAAAAAAA");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 0, 10, 9, 19, 10, "10=", 10);
    }

    #[rstest]
    fn test_local_prefer_last_jump_to_prefix_clip() {
        // NB: 10 Cs match at the end of x and start of y, so enough matches to jump after
        // matching the 10 As at the start of x and end of x.
        //    [10 .... 19] [0 ...... 9]
        // x: [CCCCCCCCCC] [AAAAAAAAAA]
        // y:  CCCCCCCCCC==>AAAAAAAAAA
        //    [0 ...... 9] [10 .... 19]
        let x = s("AAAAAAAAAACCCCCCCCCC");
        let y = s("CCCCCCCCCCAAAAAAAAAA");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 10, 10, 0, 20, 20 - 10, "10=20j10=", 20);
    }

    #[rstest]
    fn test_local_double_jump_with_trailing_y() {
        //    [0 ...... 9] [20 .... 29] [10 .... 19]
        // x: [AAAAAAAAAA] [CCCCCCCCCC] [GGGGGGGGGG] ----------
        // y:  AAAAAAAAAA==>CCCCCCCCCC==>GGGGGGGGGG  TTTTTTTTTT
        //    [0 ...... 9] [10 .... 19] [20 .... 29]
        let x = s("AAAAAAAAAAGGGGGGGGGGCCCCCCCCCC");
        let y = s("AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(
            &alignment,
            0,
            20,
            0,
            30,
            30 - 10 - 10,
            "10=10J10=20j10=",
            30,
        );
    }

    #[rstest]
    fn test_local_double_jump_with_leading_y() {
        //               [0 ...... 9] [20 .... 29] [10 .... 19]
        // x:  ----------[AAAAAAAAAA] [CCCCCCCCCC] [GGGGGGGGGG]
        // y: TTTTTTTTTT  AAAAAAAAAA==>CCCCCCCCCC==>GGGGGGGGGG
        //               [10 .... 19] [20 .... 29] [30 .... 39]
        let x = s("          AAAAAAAAAAGGGGGGGGGGCCCCCCCCCC");
        let y = s("TTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGG");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(
            &alignment,
            0,
            20,
            10,
            40,
            30 - 10 - 10,
            "10=10J10=20j10=",
            30,
        );
    }

    #[rstest]
    fn test_global_start_with_jump() {
        let x = s("TTTTTTTTTTAAAAAAAAAA");
        let y = s("          AAAAAAAAAA");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 20, 0, 10, 10 - 10, "10J10=", 10);
    }

    #[rstest]
    fn test_global_end_with_jump() {
        let x = s("AAAAAAAAAATTTTTTTTTT");
        let y = s("AAAAAAAAAA");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 20, 0, 10, 10 - 10, "10=10J", 10);
    }

    #[rstest]
    fn test_global_start_and_end_with_jump() {
        let x = s("TTTTTTTTTTAAAAAAAAAATTTTTTTTTT");
        let y = s("          AAAAAAAAAA");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 30, 0, 10, 10 - 10 - 10, "10J10=10J", 10);
    }

    #[rstest]
    fn test_local_jump_with_x_and_y() {
        let x = s("AGCT");
        let y = s("ACGT");
        // disallows mismatches and gaps, but allows jumps
        let match_fn = MatchParams::new(1, -100_000);
        let mut aligner: SingleContigAligner<MatchParams> =
            SingleContigAligner::new(-100_000, -100_000, -1, match_fn);
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 0, 4, 0, 4, 4 - 1 - 1 - 1, "1=1J1=2j1=1J1=", 4);
    }

    #[rstest]
    fn test_local_jump_with_x_and_y_suffix_clips_small() {
        let x = s("AAGGCCT");
        let y = s("AACCGGT");
        // disallows mismatches and gaps, but allows jumps
        let match_fn = MatchParams::new(1, -100_000);
        let mut aligner = SingleContigAligner::new(-100_000, -100_000, -2, match_fn);
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 0, 4, 0, 6, 6 - 2 - 2, "2=2J2=4j2=", 6);
    }

    #[rstest]
    fn test_local_jump_with_x_and_y_suffix_clips() {
        //    [0 ...... 9] [20 .... 29] [10 .... 19]
        // x: [AAAAAAAAAA] [CCCCCCCCCC] [GGGGGGGGGG]          TTTTTTTTT
        // y:  AAAAAAAAAA==>CCCCCCCCCC==>GGGGGGGGGG TTTTTTTTT
        //    [0 ...... 9] [10 .... 19] [20 .... 29]
        let x = s("AAAAAAAAAAGGGGGGGGGGCCCCCCCCCCTTTTTTTTT");
        let y = s("AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTT");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(
            &alignment,
            0,
            20,
            0,
            30,
            30 - 10 - 10,
            "10=10J10=20j10=",
            30,
        );
    }

    #[rstest]
    fn test_local_jump_with_x_and_y_prefix_clips_small() {
        let x = s("AGGCCTT");
        let y = s("ACCGGTT");
        // disallows mismatches and gaps, but allows jumps
        let match_fn = MatchParams::new(1, -100_000);
        let mut aligner = SingleContigAligner::new(-100_000, -100_000, -2, match_fn);
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 3, 7, 1, 7, 6 - 2 - 2, "2=4j2=2J2=", 6);
    }

    #[rstest]
    fn test_local_jump_with_x_and_y_prefix_clips() {
        //              [20 .....29] [10 .... 19] [30 .... 39]
        // x: TTTTTTTTT [GGGGGGGGGG] [CCCCCCCCCC] [AAAAAAAAAA]
        // y: TTTTTTTTT  GGGGGGGGGG==>CCCCCCCCCC==>AAAAAAAAAA
        //              [10 .... 19] [20 .... 29] [30 .... 39]
        let x = s("TTTTTTTTTCCCCCCCCCCGGGGGGGGGGAAAAAAAAAA");
        let y = s("TTTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAA");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        let alignment = aligner.local(&x, &y);
        assert_alignment(
            &alignment,
            19,
            39,
            9,
            39,
            30 - 10 - 10,
            "10=20j10=10J10=",
            30,
        );
    }

    #[rstest]
    fn test_local_jump() {
        //              [20 .....29] [10 .... 19] [30 .... 39]
        // x: TTTTTTTTT [GGGGGGGGGG] [CCCCCCCCCC] [AAAAAAAAAA]
        // y: TTTTTTTTT  GGGGGGGGGG==>CCCCCCCCCC==>AAAAAAAAAA
        //              [10 .... 19] [20 .... 29] [30 .... 39]
        let x = s("TTTTTTTTTCCCCCCCCCCGGGGGGGGGGAAAAAAAAAA");
        let y = s("TTTTTTTTTGGGGGGGGGGCCCCCCCCCCAAAAAAAAAA");
        let mut aligner: SingleContigAligner<MatchParams> = SingleContigAligner::default();
        aligner.scoring = aligner.scoring.set_jump_score(-10);
        let alignment = aligner.local(&x, &y);
        assert_alignment(
            &alignment,
            19,
            39,
            9,
            39,
            30 - 10 - 10,
            "10=20j10=10J10=",
            30,
        );
    }

    #[rstest]
    fn test_global_short_jumps() {
        let x = s("AAGGCCTT");
        let y = s("AACCGGTT");
        let match_fn = MatchParams::new(1, -100_000);
        let mut aligner = SingleContigAligner::new(-100_000, -100_000, -1, match_fn);
        let alignment = aligner.global(&x, &y);
        assert_alignment(&alignment, 0, 8, 0, 8, 8 - 1 - 1 - 1, "2=2J2=4j2=2J2=", 8);
    }

    #[rstest]
    fn test_local_circular_jump() {
        let x = s("AACCGGTT");
        let y = s("TTAA");
        let match_fn = MatchParams::new(1, -100_000);
        let mut aligner = SingleContigAligner::new(-100_000, -100_000, -1, match_fn);
        aligner.set_circular(true);
        let alignment = aligner.local(&x, &y);
        assert_alignment(&alignment, 6, 2, 0, 4, 4, "2=8j2=", 4);
    }

    #[rstest]
    fn test_targetlocal_circular_jump() {
        let x = s("GGTTAACC");
        let y = s("AACCGGTT");
        let match_fn = MatchParams::new(1, -100_000);
        let mut aligner = SingleContigAligner::new(-100_000, -100_000, -1, match_fn);
        aligner.set_circular(true);
        let alignment = aligner.targetlocal(&x, &y);
        assert_alignment(&alignment, 4, 4, 0, 8, 8, "4=8j4=", 8);
    }
}
