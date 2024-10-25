//! Compute all combinations of values in a given list
//! Credit: https://github.com/meltinglava/uniquecombinations/
#[cfg(not(feature = "std"))]
use ark_std::vec::Vec;
/// Compute all combinations of values in a given list.
pub(crate) struct Combinations<T>
where
    T: Ord + Clone,
{
    original: Vec<T>,
    possition: Vec<usize>,
    len: usize,
    started: bool,
}

impl<T> Combinations<T>
where
    T: Ord + Clone,
{
    /// Initialize the permutations.
    pub(crate) fn new(mut original: Vec<T>, len: usize) -> Self {
        if original.len() > len && len >= 1 {
            original.sort_unstable();
            Self {
                original,
                possition: (0..len).collect(),
                len,
                started: false,
            }
        } else {
            panic!("the length has to be smaller then the datasets len");
        }
    }

    #[inline]
    fn insert(&self, col: &mut Vec<T>) {
        col.clear();
        self.possition
            .iter()
            .enumerate()
            .for_each(|(p, n)| col.insert(p, self.original[*n].clone()))
    }

    /// Clear the contents of the comb vector and insert the next combination.
    fn next_combination(&mut self, mut comb: &mut Vec<T>) -> bool {
        if !self.started {
            // first pass throught
            self.started = true;
            self.insert(&mut comb);
            true
        } else {
            let org_len = self.original.len();
            // check if we cant bump the back number
            if self.original[self.possition[self.len - 1]] == self.original[org_len - 1] {
                // locate the number closest behind that needs to be bumped
                for i in 2..=self.len {
                    if self.original[self.possition[self.len - i]] < self.original[org_len - i] {
                        //find the value of the
                        let lastpos = self.possition[self.len - i];
                        let val = &self.original[lastpos];
                        for j in lastpos + 1..org_len {
                            if *val < self.original[j] {
                                for k in 0..i {
                                    self.possition[self.len - i + k] = j + k;
                                }
                                self.insert(&mut comb);
                                return true;
                            }
                        }
                    }
                }
                false
            } else {
                let mut i = self.possition[self.len - 1];
                let current = &self.original[i];
                let mut next = current;
                while current == next {
                    i += 1;
                    next = &self.original[i];
                }
                self.possition[self.len - 1] = i;
                self.insert(&mut comb);
                true
            }
        }
    }
}

impl<T> Iterator for Combinations<T>
where
    T: Ord + Clone,
{
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut vals = Vec::with_capacity(self.len);
        if self.next_combination(&mut vals) {
            Some(vals)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equals() {
        assert!(Combinations::new(vec![2, 2, 2], 2).next().unwrap() == vec![2, 2])
    }

    #[test]
    fn t_123() {
        assert!(
            dbg!(Combinations::new(vec![1, 2, 3], 2)
                .take(10)
                .collect::<Vec<_>>())
                == vec![vec![1, 2], vec![1, 3], vec![2, 3]]
        )
    }

    #[test]
    fn complicated() {
        let actual: Vec<_> = Combinations::new(vec![1, 2, 2, 3, 4], 3).collect();
        let expected = vec![
            vec![1, 2, 2],
            vec![1, 2, 3],
            vec![1, 2, 4],
            vec![1, 3, 4],
            vec![2, 2, 3],
            vec![2, 2, 4],
            vec![2, 3, 4],
        ];
        assert!(actual == expected)
    }
}
