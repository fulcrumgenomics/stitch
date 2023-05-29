use serde::Deserialize;
use serde::Serialize;

use super::Cell;
use super::TracebackCell;
use super::TB_START;

/// Internal traceback.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Traceback {
    rows: usize,
    cols: usize,
    matrix: Vec<Cell>,
}

impl Traceback {
    pub fn with_capacity(m: usize, n: usize) -> Self {
        let rows = m + 1;
        let cols = n + 1;
        Traceback {
            rows,
            cols,
            matrix: Vec::with_capacity(rows * cols),
        }
    }

    pub fn init(&mut self, m: usize, n: usize) {
        self.matrix.clear();
        let mut start = crate::alignment::traceback::default();
        start.set_all(TB_START, 0);
        start.set_s_all(TB_START, 0, 0, false);
        // set every cell to start
        self.resize(m, n, start);
    }

    #[inline(always)]
    pub fn set(&mut self, i: usize, j: usize, v: Cell) {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        self.matrix[i * self.cols + j] = v;
    }

    #[inline(always)]
    pub fn get(&self, i: usize, j: usize) -> &Cell {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        &self.matrix[i * self.cols + j]
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut Cell {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        &mut self.matrix[i * self.cols + j]
    }

    pub fn resize(&mut self, m: usize, n: usize, v: Cell) {
        self.rows = m + 1;
        self.cols = n + 1;
        self.matrix.resize(self.rows * self.cols, v);
    }
}
