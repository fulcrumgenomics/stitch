// Original Code was copied from:
//     https://github.com/rust-bio/rust-bio/blob/master/src/alignment/pairwise/mod.rs
// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp::max;
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
            self.from[i] = from;
        } else {
            self.from[i] = from + 0xFFFF_FFFF;
        }
    }

    #[inline(always)]
    pub fn get(&self, i: usize) -> (i32, u32, bool) {
        (
            self.score[i],
            self.from[i] << 1 >> 1,
            (self.from[i] >> 31) == 0,
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
        // NB: do not jump from the first/last `i`
        self.score[m] = MIN_SCORE;
        self.from[m] = m as u32;
        for i in (1..m).rev() {
            if self.score[i + 1] >= S[prev][i] {
                self.score[i] = self.score[i + 1];
                self.from[i] = self.from[i + 1];
            } else {
                self.score[i] = S[prev][i];
                self.from[i] = i as u32;
            };
        }
        self.score[0] = MIN_SCORE;

        let mut earlier_jump_score = MIN_SCORE;
        let mut earlier_jump_length = 0u32;
        self.score[0] += scoring.jump_score;
        for i in 1..=m {
            if 2 <= i {
                let score = S[prev][i - 2] + scoring.jump_score;
                if score >= earlier_jump_score {
                    earlier_jump_score = score;
                    earlier_jump_length = (i - 2) as u32;
                }
            }

            let later_jump_score = self.score[i] + scoring.jump_score;
            let diagonal_score = S[prev][i - 1];
            self.score[i] = max(earlier_jump_score, max(diagonal_score, later_jump_score));
            if diagonal_score == self.score[i] {
                self.from[i] = (i - 1) as u32;
            } else if earlier_jump_score == self.score[i] {
                self.from[i] = earlier_jump_length;
            }
        }
    }
}
