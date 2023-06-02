use bio::alignment::pairwise::MatchFunc;
use serde::Serialize;

use crate::align::aligners::constants::MIN_SCORE;

/// Details of scoring are encapsulated in this structure.
///
/// An [affine gap score model](https://en.wikipedia.org/wiki/Gap_penalty#Affine)
/// is used so that the gap score for a length `k` is:
/// `GapScore(k) = gap_open + gap_extend * k`
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct Scoring<F: MatchFunc> {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub jump_score_same_contig_and_strand: i32,
    pub jump_score_same_contig_opposite_strand: i32,
    pub jump_score_inter_contig: i32,
    pub match_fn: F,
    pub match_scores: Option<(i32, i32)>,
    pub xclip_prefix: i32,
    pub xclip_suffix: i32,
    pub yclip_prefix: i32,
    pub yclip_suffix: i32,
}

impl<F: MatchFunc> Scoring<F> {
    /// Create new Scoring instance with given gap open, gap extend penalties
    /// and the score function. The clip penalties are set to [`MIN_SCORE`](constant.MIN_SCORE.html) by default
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `jump_score` - the score for jumping in the query (should not be positive)
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    pub fn with_jump_score(gap_open: i32, gap_extend: i32, jump_score: i32, match_fn: F) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");
        assert!(jump_score <= 0, "jump_score can't be positive");

        Self::with_jump_scores(
            gap_open, gap_extend, jump_score, jump_score, jump_score, match_fn,
        )
    }

    /// Create new Scoring instance with given gap open, gap extend penalties
    /// and the score function. The clip penalties are set to [`MIN_SCORE`](constant.MIN_SCORE.html) by default
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `jump_score_same_contig_and_strand` - the score for jumping to the same contig and strand in the query (should not be positive)
    /// * `jump_score_same_contig_opposite_strand` - the score for jumping to the same contig and opposite strand in the query (should not be positive)
    /// * `jump_score_inter_contig` - the score for jumping to a different contig in the query (should not be positive)
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    pub fn with_jump_scores(
        gap_open: i32,
        gap_extend: i32,
        jump_score_same_contig_and_strand: i32,
        jump_score_same_contig_opposite_strand: i32,
        jump_score_inter_contig: i32,
        match_fn: F,
    ) -> Self {
        assert!(
            jump_score_same_contig_and_strand <= 0,
            "jump_score_same_contig_and_strand can't be positive"
        );
        assert!(
            jump_score_same_contig_opposite_strand <= 0,
            "jump_score_same_contig_opposite_strand can't be positive"
        );
        assert!(
            jump_score_inter_contig <= 0,
            "jump_score_inter_contig can't be positive"
        );

        Self {
            gap_open,
            gap_extend,
            jump_score_same_contig_and_strand,
            jump_score_same_contig_opposite_strand,
            jump_score_inter_contig,
            match_fn,
            match_scores: None,
            xclip_prefix: MIN_SCORE,
            xclip_suffix: MIN_SCORE,
            yclip_prefix: MIN_SCORE,
            yclip_suffix: MIN_SCORE,
        }
    }

    /// Sets the jump scores to the given value
    ///
    /// # Arguments
    ///
    /// * `jump_score` - Jump score
    #[allow(dead_code)]
    pub fn set_jump_score(mut self, jump_score: i32) -> Self {
        self.jump_score_same_contig_and_strand = jump_score;
        self.jump_score_same_contig_opposite_strand = jump_score;
        self.jump_score_inter_contig = jump_score;
        self
    }

    /// Sets the jump scores to the given values
    #[allow(dead_code)]
    pub fn set_jump_scores(
        mut self,
        jump_score_same_contig_and_strand: i32,
        jump_score_same_contig_opposite_strand: i32,
        jump_score_inter_contig: i32,
    ) -> Self {
        self.jump_score_same_contig_and_strand = jump_score_same_contig_and_strand;
        self.jump_score_same_contig_opposite_strand = jump_score_same_contig_opposite_strand;
        self.jump_score_inter_contig = jump_score_inter_contig;
        self
    }

    /// Sets the prefix and suffix clipping penalties for x to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Clipping penalty for x (both prefix and suffix, should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).xclip(-5);
    /// assert!(scoring.xclip_prefix == -5);
    /// assert!(scoring.yclip_prefix == MIN_SCORE);
    /// assert!(scoring.xclip_suffix == -5);
    /// assert!(scoring.yclip_suffix == MIN_SCORE);
    /// ```
    #[allow(dead_code)]
    pub fn xclip(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.xclip_prefix = penalty;
        self.xclip_suffix = penalty;
        self
    }

    /// Sets the prefix clipping penalty for x to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Prefix clipping penalty for x (should not be positive)
    ///
    /// # Example
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).xclip_prefix(-5);
    /// assert!(scoring.xclip_prefix == -5);
    /// assert!(scoring.yclip_prefix == MIN_SCORE);
    /// assert!(scoring.xclip_suffix == MIN_SCORE);
    /// assert!(scoring.yclip_suffix == MIN_SCORE);
    /// ```
    #[allow(dead_code)]
    pub fn xclip_prefix(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.xclip_prefix = penalty;
        self
    }

    /// Sets the suffix clipping penalty for x to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Suffix clipping penalty for x (should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).xclip_suffix(-5);
    /// assert!(scoring.xclip_prefix == MIN_SCORE);
    /// assert!(scoring.yclip_prefix == MIN_SCORE);
    /// assert!(scoring.xclip_suffix == -5);
    /// assert!(scoring.yclip_suffix == MIN_SCORE);
    /// ```
    #[allow(dead_code)]
    pub fn xclip_suffix(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.xclip_suffix = penalty;
        self
    }

    /// Sets the prefix and suffix clipping penalties for y to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Clipping penalty for y (both prefix and suffix, should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).yclip(-5);
    /// assert!(scoring.xclip_prefix == MIN_SCORE);
    /// assert!(scoring.yclip_prefix == -5);
    /// assert!(scoring.xclip_suffix == MIN_SCORE);
    /// assert!(scoring.yclip_suffix == -5);
    /// ```
    #[allow(dead_code)]
    pub fn yclip(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.yclip_prefix = penalty;
        self.yclip_suffix = penalty;
        self
    }

    /// Sets the prefix clipping penalty for y to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Prefix clipping penalty for y (should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).yclip_prefix(-5);
    /// assert!(scoring.xclip_prefix == MIN_SCORE);
    /// assert!(scoring.yclip_prefix == -5);
    /// assert!(scoring.xclip_suffix == MIN_SCORE);
    /// assert!(scoring.yclip_suffix == MIN_SCORE);
    /// ```
    #[allow(dead_code)]
    pub fn yclip_prefix(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.yclip_prefix = penalty;
        self
    }

    /// Sets the suffix clipping penalty for y to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Suffix clipping penalty for y (should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).yclip_suffix(-5);
    /// assert!(scoring.xclip_prefix == MIN_SCORE);
    /// assert!(scoring.yclip_prefix == MIN_SCORE);
    /// assert!(scoring.xclip_suffix == MIN_SCORE);
    /// assert!(scoring.yclip_suffix == -5);
    /// ```
    #[allow(dead_code)]
    pub fn yclip_suffix(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.yclip_suffix = penalty;
        self
    }
}
