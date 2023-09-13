/// A `HashMap` where the keys are `usize`s and values can be anything.
///
/// Uses a `Vec` internally to store the values, up to the maximum key size seen so far.  This
/// provides performance for get, put, and contain operations.
pub(crate) struct IndexMap<T> {
    data: Vec<Option<T>>,
    total_added: usize,
}

#[allow(dead_code)]
impl<T: Clone> IndexMap<T> {
    pub fn new(max_index: usize) -> Self {
        Self {
            data: (0..=max_index)
                .map(|_| None::<T>)
                .collect::<Vec<Option<T>>>(),
            total_added: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.total_added
    }

    pub fn is_empty(&self) -> bool {
        self.total_added == 0
    }

    pub fn keys(&self) -> impl Iterator<Item = usize> + '_ {
        self.data
            .iter()
            .enumerate()
            .filter_map(|(i, v)| if v.is_none() { None } else { Some(i) })
    }

    pub fn values(&self) -> impl Iterator<Item = &T> + '_ {
        self.data.iter().flatten()
    }

    pub fn contains(&self, index: usize) -> bool {
        self.data[index].is_some()
    }

    #[inline(always)]
    pub fn contains_u32(&self, index: u32) -> bool {
        self.contains(index as usize)
    }

    pub fn get(&self, index: usize) -> Option<T> {
        assert!(index < self.data.len());
        self.data[index].clone()
    }

    #[inline(always)]
    pub fn get_u32(&self, index: u32) -> Option<T> {
        self.get(index as usize)
    }

    pub fn put(&mut self, index: usize, value: T) {
        if self.data.len() <= index {
            self.reserve(index + 1);
        }
        if self.data[index].is_none() {
            self.total_added += 1;
        }
        self.data[index] = Some(value);
    }

    #[inline(always)]
    pub fn put_u32(&mut self, index: u32, value: T) {
        self.put(index as usize, value);
    }

    pub fn capacity(&self) -> usize {
        self.data.len()
    }

    pub fn reserve(&mut self, capacity: usize) {
        while self.capacity() < capacity {
            self.data.push(None);
        }
    }
}

// Tests
#[cfg(test)]
pub mod tests {
    use itertools::Itertools;
    use rstest::rstest;

    use crate::util::index_map::IndexMap;

    #[rstest]
    fn test_new_empty() {
        let imap = IndexMap::<usize>::new(12);
        assert_eq!(imap.len(), 0);
        assert!(imap.is_empty());
        assert!(imap.keys().collect_vec().is_empty());
        assert!(imap.values().collect_vec().is_empty());
        assert_eq!(imap.capacity(), 13);
    }

    #[rstest]
    fn test_put_and_get() {
        let mut imap = IndexMap::<usize>::new(4);

        imap.put(1, 2);
        assert_eq!(imap.len(), 1);
        assert!(!imap.is_empty());
        assert_eq!(imap.keys().collect_vec(), vec![1]);
        assert_eq!(imap.values().collect_vec(), vec![&2]);
        assert_eq!(imap.capacity(), 5);
        assert_eq!(imap.get(0), None);
        assert_eq!(imap.get(1), Some(2));
        assert_eq!(imap.get(2), None);
        assert_eq!(imap.get(3), None);

        imap.put(10, 4);
        assert_eq!(imap.len(), 2);
        assert!(!imap.is_empty());
        assert_eq!(imap.keys().collect_vec(), vec![1, 10]);
        assert_eq!(imap.values().collect_vec(), vec![&2, &4]);
        assert_eq!(imap.capacity(), 11);
        assert_eq!(imap.get(0), None);
        assert_eq!(imap.get(1), Some(2));
        assert_eq!(imap.get(2), None);
        assert_eq!(imap.get(3), None);
        assert_eq!(imap.get(10), Some(4));

        imap.reserve(11);
        assert_eq!(imap.capacity(), 11);

        imap.reserve(12);
        assert_eq!(imap.capacity(), 12);
        assert_eq!(imap.get(11), None);

        imap.put(11, 124);
        assert_eq!(imap.get(11), Some(124));
    }
}
