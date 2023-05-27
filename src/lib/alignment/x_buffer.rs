// Original Code was copied from:
//     https://github.com/rust-bio/rust-bio/blob/master/src/alignment/pairwise/mod.rs
// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::i32;
use std::iter::repeat;

use bio::alignment::pairwise::MatchFunc;
use serde::Deserialize;
use serde::Serialize;

use crate::alignment::scoring::Scoring;

use super::constants::MIN_SCORE;

#[allow(non_snake_case)]
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct XBuffer {
    score: Vec<i32>,
    from: Vec<u32>,
}

impl XBuffer {
    pub fn new(m: usize) -> Self {
        let mut score_buffer = Vec::with_capacity(m + 1);
        let mut from_buffer = Vec::with_capacity(m + 1);
        score_buffer.extend(repeat(MIN_SCORE).take(m + 1));
        from_buffer.extend(repeat(0).take(m + 1));
        XBuffer {
            score: score_buffer,
            from: from_buffer,
        }
    }

    #[inline(always)]
    pub fn set(&mut self, i: usize, score: i32, from: u32, flip_strand: bool) {
        self.score[i] = score;
        if flip_strand {
            self.from[i] = from + (1 << 31);
        } else {
            self.from[i] = from;
        }
    }

    #[inline(always)]
    pub fn get(&self, i: usize) -> (i32, u32, bool) {
        (
            self.score[i],
            self.from[i] << 1 >> 1,
            (self.from[i] >> 31) == 1,
        )
    }

    #[allow(non_snake_case)]
    pub fn fill<F: MatchFunc>(
        &mut self,
        m: usize,
        prev: usize,
        S: &[Vec<i32>; 2],
        scoring: &Scoring<F>,
    ) {
        // The current base in `y` is fixed (i.e. index `j`), but we can move from anywhere in `x`.
        // For a given `i`, we want to find jump from `k` that yields the maximum score.
        // We can compute this in two passes:
        // 1. suffix_score[k] = max(suffix_score[k + 1], S[prev][k] + jump_penalty)
        // 2. prefix_score[k] = max(prefix_score[k - 1], S[prev][k] + jump_penalty)
        // Then score[k] = max(prefix_score[k], suffix_score[k]).  This can be computer at the
        // same time as prefix_score[k].

        // The first pass to compute for each k:
        //   suffix_score[k] = max(suffix_score[k + 1], S[prev][k] + jump_penalty)
        // base case
        self.score[m] = S[prev][m] + scoring.jump_score;
        self.from[m] = m as u32;
        // iterate over k in descending order
        for k in (0..m).rev() {
            let k_score = S[prev][k] + scoring.jump_score;
            if self.score[k + 1] >= k_score {
                self.score[k] = self.score[k + 1];
                self.from[k] = self.from[k + 1];
            } else {
                self.score[k] = k_score;
                self.from[k] = k as u32;
            };
        }

        // The second pass to compute for each k:
        //   prefix_score[k] = max(prefix_score[k - 1], S[prev][k] + jump_penalty)
        // and then
        //   score[k] = max(prefix_score[k], suffix_score[k])
        // base case
        let mut prev_prefix_score = S[prev][0] + scoring.jump_score;
        let mut prev_prefix_from = 0u32;
        // iterate over k in ascending order
        for k in 1..=m {
            // compute prefix_score[k]
            let cur_prefix_score = S[prev][k] + scoring.jump_score;
            if cur_prefix_score >= prev_prefix_score {
                prev_prefix_score = cur_prefix_score;
                prev_prefix_from = k as u32;
            }
            // score[k] = max(prefix_score[k], suffix_score[k])
            if prev_prefix_score >= self.score[k] {
                self.score[k] = prev_prefix_score;
                self.from[k] = prev_prefix_from;
            }
        }
    }
}
