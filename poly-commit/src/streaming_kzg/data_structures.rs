use ark_ff::Field;
use ark_std::borrow::Borrow;
use ark_std::vec::Vec;

use crate::streaming_kzg::ceil_div;
use ark_std::iterable::Iterable;

/// A `Streamer` folding a vector of coefficients
/// with the given challenges, and producing a stream of items
/// `(i, v)` where `i` indicates the depth, and `v` is the next coefficient.
/// The stream can produce all foldings in the tree with a single pass over the initial stream.
#[derive(Clone, Copy)]
pub struct FoldedPolynomialTree<'a, F, S> {
    challenges: &'a [F],
    coefficients: &'a S,
}

impl<'a, F, S> FoldedPolynomialTree<'a, F, S>
where
    S: Iterable,
    F: Field,
    S::Item: Borrow<F>,
{
    /// Initialize a new polynomial tree.
    pub fn new(coefficients: &'a S, challenges: &'a [F]) -> Self {
        Self {
            coefficients,
            challenges,
        }
    }

    /// Outputs the depth of the polynomial tree.
    #[inline]
    pub fn depth(&self) -> usize {
        self.challenges.len()
    }
}

impl<'a, F, S> Iterable for FoldedPolynomialTree<'a, F, S>
where
    S: Iterable,
    F: Field,
    S::Item: Borrow<F>,
{
    type Item = (usize, F);

    type Iter = FoldedPolynomialTreeIter<'a, F, S::Iter>;

    fn iter(&self) -> Self::Iter {
        FoldedPolynomialTreeIter::new(
            self.coefficients.iter(),
            self.coefficients.len(),
            self.challenges,
        )
    }

    fn len(&self) -> usize {
        self.coefficients.len()
    }
}

/// Iterator of the polynomial tree.
pub struct FoldedPolynomialTreeIter<'a, F, I> {
    challenges: &'a [F],
    iterator: I,
    stack: Vec<(usize, F)>,
}

fn init_stack<F: Field>(n: usize, challenges_len: usize) -> Vec<(usize, F)> {
    let mut stack = Vec::with_capacity(challenges_len);

    // generally we expect the size to be a power of two.
    // If not, we are going to fill the stack as if the array was padded to zero up to the expected size.
    let chunk_size = 1 << challenges_len;
    if n % chunk_size != 0 {
        let mut delta = chunk_size - n % chunk_size;
        for i in (0..challenges_len).rev() {
            if delta >= 1 << i {
                stack.push((i, F::zero()));
                delta -= 1 << i
            }
        }
    }
    stack
}

impl<'a, F, I> FoldedPolynomialTreeIter<'a, F, I>
where
    F: Field,
    I: Iterator,
    I::Item: Borrow<F>,
{
    fn new(iterator: I, n: usize, challenges: &'a [F]) -> Self {
        let stack = init_stack(n, challenges.len());

        Self {
            challenges,
            iterator,
            stack,
        }
    }
}

impl<'a, F, I> Iterator for FoldedPolynomialTreeIter<'a, F, I>
where
    F: Field,
    I: Iterator,
    I::Item: Borrow<F>,
{
    type Item = (usize, F);

    fn next(&mut self) -> Option<<Self as Iterator>::Item> {
        let len = self.stack.len();
        let item = if len > 1 && self.stack[len - 1].0 == self.stack[len - 2].0 {
            // pop the last two elements from the stack.
            // we could also use .pop() twice but truncate is slightly faster.
            let (_level, lhs) = self.stack[len - 1];
            let (level, rhs) = self.stack[len - 2];
            self.stack.truncate(len - 2);
            // fold them producing the coefficient and the level `level+1`
            let folded_coefficient = rhs * self.challenges[level] + lhs;
            (level + 1, folded_coefficient)
        } else {
            (0, *self.iterator.next()?.borrow())
        };

        // do not add to the stack the coefficient of the max-depth folded polynomial.
        if item.0 != self.challenges.len() {
            self.stack.push(item)
        }

        // Skip the base polynomial, recursively calling itself to access the next level
        if item.0 == 0 {
            self.next()
        } else {
            Some(item)
        }
    }
}

/// Stream implementation of foleded polynomial.
#[derive(Clone, Copy)]
pub struct FoldedPolynomialStream<'a, F, S>(FoldedPolynomialTree<'a, F, S>, usize);
/// Iterator implementation of foleded polynomial.
pub struct FoldedPolynomialStreamIter<'a, F, I> {
    challenges: &'a [F],
    iterator: I,
    stack: Vec<(usize, F)>,
}

impl<'a, F, S> FoldedPolynomialStream<'a, F, S>
where
    S: Iterable,
    F: Field,
    S::Item: Borrow<F>,
{
    /// Initialize a new folded polynomial stream.
    pub fn new(coefficients: &'a S, challenges: &'a [F]) -> Self {
        let tree = FoldedPolynomialTree::new(coefficients, challenges);
        let len = challenges.len();
        Self(tree, len)
    }
}

impl<'a, F, S> Iterable for FoldedPolynomialStream<'a, F, S>
where
    S: Iterable,
    F: Field,
    S::Item: Borrow<F>,
{
    type Item = F;
    type Iter = FoldedPolynomialStreamIter<'a, F, S::Iter>;

    fn iter(&self) -> Self::Iter {
        let iterator = self.0.coefficients.iter();
        let challenges = self.0.challenges;
        let stack = init_stack(self.0.coefficients.len(), challenges.len());
        FoldedPolynomialStreamIter {
            iterator,
            challenges,
            stack,
        }
    }

    fn len(&self) -> usize {
        ceil_div(self.0.len(), 1 << self.0.challenges.len())
    }
}

impl<'a, F, I> Iterator for FoldedPolynomialStreamIter<'a, F, I>
where
    F: Field,
    I: Iterator,
    I::Item: Borrow<F>,
{
    type Item = F;

    fn next(&mut self) -> Option<Self::Item> {
        let target_level = self.challenges.len();
        loop {
            let len = self.stack.len();
            let (level, element) = if len > 1 && self.stack[len - 1].0 == self.stack[len - 2].0 {
                let (_level, lhs) = self.stack[len - 1];
                let (level, rhs) = self.stack[len - 2];
                self.stack.truncate(len - 2);

                let folded_coefficient = rhs * self.challenges[level] + lhs;
                (level + 1, folded_coefficient)
            } else if target_level > 0 && (len == 0 || (len > 0 && self.stack[len - 1].0 != 0)) {
                // If the target level is strictly positive, there's no need to put elements of level zero in the stream.
                // We can immediately read 2 elements from the stream and push an element of the form (1, folded_coefficient).
                // Nota bene: this branch is not needed, but brings in a decent speed-up for the resulting implementation.
                let rhs = self.iterator.next()?;
                let lhs = self.iterator.next()?;

                let folded_coefficient = self.challenges[0] * rhs.borrow() + lhs.borrow();
                (1, folded_coefficient)
            } else {
                (0, *self.iterator.next()?.borrow())
            };

            // do not add to the stack the coefficient of the folded polynomial, but instead return it.
            if level != target_level {
                self.stack.push((level, element))
            } else {
                return Some(element);
            }
        }
    }
}

#[test]
fn test_folded_polynomial() {
    use ark_bls12_381::Fr as F;
    use ark_ff::One;

    let two = F::one() + F::one();

    let coefficients = vec![F::one(), two, F::one(), F::one()];
    let challenges = vec![F::one(), two];
    let coefficients_stream = coefficients.as_slice();
    let foldstream = FoldedPolynomialTree::new(&coefficients_stream, challenges.as_slice());
    let fold_stream = FoldedPolynomialStream(foldstream, 2);
    assert_eq!(fold_stream.len(), 1);
    assert_eq!(
        fold_stream.iter().next(),
        Some(two + two * (F::one() + two))
    );

    let one = F::one();
    let coefficients = vec![one; 12];
    let challenges = vec![F::one(); 4];
    let coefficients_stream = coefficients.as_slice();
    let foldstream = FoldedPolynomialTree::new(&coefficients_stream, challenges.as_slice());
    let fold_stream = FoldedPolynomialStream(foldstream, 4).iter();
    assert_eq!(fold_stream.last(), Some(coefficients.iter().sum()));
}

#[test]
fn test_folded_polynomial_tree() {
    use ark_bls12_381::Fr as F;
    use ark_ff::One;

    let two = F::one() + F::one();

    let coefficients = vec![F::one(), two, F::one(), F::one()];
    let challenges = vec![F::one(), two];
    let coefficients_stream = coefficients.as_slice();
    let fold_streamer = FoldedPolynomialTree::new(&coefficients_stream, challenges.as_slice());
    let mut fold_iter = fold_streamer.iter();
    // assert_eq!(fold_stream.next(), Some((0, F::one())));
    // assert_eq!(fold_stream.next(), Some((0, two)));
    assert_eq!(fold_iter.next(), Some((1, F::one() + two)));
    // assert_eq!(fold_stream.next(), Some((0, F::one())));
    // assert_eq!(fold_stream.next(), Some((0, F::one())));
    assert_eq!(fold_iter.next(), Some((1, F::one() + F::one())));
    assert_eq!(fold_iter.next(), Some((2, two + two * (F::one() + two))));

    let one = F::one();
    let coefficients = vec![one; 12];
    let challenges = vec![F::one(); 4];
    let coefficients_stream = coefficients.as_slice();
    let fold_streamer = FoldedPolynomialTree::new(&coefficients_stream, challenges.as_slice());
    let fold_init = fold_streamer.iter();
    let mut fold_iter = fold_init.skip(5);
    assert_eq!(fold_iter.next(), Some((1, two)));
    assert_eq!(fold_iter.last(), Some((4, coefficients.iter().sum())));
}
