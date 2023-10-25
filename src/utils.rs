use core::marker::PhantomData;

#[cfg(not(feature = "std"))]
use num_traits::Float;

#[cfg(feature = "parallel")]
use rayon::{
    iter::{IntoParallelRefIterator, ParallelIterator},
    prelude::IndexedParallelIterator,
};

use ark_ff::{Field, PrimeField};
use ark_serialize::CanonicalSerialize;
use ark_std::vec::Vec;
use merlin::Transcript;

use crate::Error;

/// Takes as input a struct, and converts them to a series of bytes. All traits
/// that implement `CanonicalSerialize` can be automatically converted to bytes
/// in this manner.
/// From jellyfish lib
#[macro_export]
macro_rules! to_bytes {
    ($x:expr) => {{
        let mut buf = ark_std::vec![];
        ark_serialize::CanonicalSerialize::serialize_compressed($x, &mut buf).map(|_| buf)
    }};
}

/// Return ceil(x / y).
pub(crate) fn ceil_div(x: usize, y: usize) -> usize {
    // XXX. warning: this expression can overflow.
    (x + y - 1) / y
}

#[derive(Debug)]
pub(crate) struct Matrix<F: Field> {
    pub(crate) n: usize,
    pub(crate) m: usize,
    entries: Vec<Vec<F>>,
}

impl<F: Field> Matrix<F> {
    /// Returns a Matrix given a list of its rows, each in turn represented as a list of field elements.
    ///
    /// # Panics
    /// Panics if the sub-lists do not all have the same length.
    pub(crate) fn new_from_rows(row_list: Vec<Vec<F>>) -> Self {
        let m = row_list[0].len();

        for row in row_list.iter().skip(1) {
            assert_eq!(
                row.len(),
                m,
                "Invalid matrix construction: not all rows have the same length"
            );
        }

        Self {
            n: row_list.len(),
            m,
            entries: row_list,
        }
    }

    /// Returns the entry in position (i, j). **Indexing starts at 0 in both coordinates**,
    /// i.e. the first element is in position (0, 0) and the last one in (n - 1, j - 1),
    /// where n and m are the number of rows and columns, respectively.
    ///
    /// Index bound checks are waived for efficiency and behaviour under invalid indexing is undefined
    #[cfg(test)]
    pub(crate) fn entry(&self, i: usize, j: usize) -> F {
        self.entries[i][j]
    }

    /// Returns the product v * self, where v is interpreted as a row vector. In other words,
    /// it returns a linear combination of the rows of self with coefficients given by v.
    ///
    /// Panics if the length of v is different from the number of rows of self.
    pub(crate) fn row_mul(&self, v: &[F]) -> Vec<F> {
        assert_eq!(
            v.len(),
            self.n,
            "Invalid row multiplication: vector has {} elements whereas each matrix column has {}",
            v.len(),
            self.n
        );

        (0..self.m)
            .map(|col| {
                inner_product(
                    v,
                    &(0..self.n)
                        .map(|row| self.entries[row][col])
                        .collect::<Vec<F>>(),
                )
            })
            .collect()
    }
}

#[inline]
pub(crate) fn inner_product<F: Field>(v1: &[F], v2: &[F]) -> F {
    ark_std::cfg_iter!(v1)
        .zip(v2)
        .map(|(li, ri)| *li * ri)
        .sum()
}

#[inline]
pub(crate) fn scalar_by_vector<F: Field>(s: F, v: &[F]) -> Vec<F> {
    ark_std::cfg_iter!(v).map(|x| *x * s).collect()
}

#[inline]
pub(crate) fn vector_sum<F: Field>(v1: &[F], v2: &[F]) -> Vec<F> {
    ark_std::cfg_iter!(v1)
        .zip(v2)
        .map(|(li, ri)| *li + ri)
        .collect()
}

/// The following struct is taken from jellyfish repository. Once they change
/// their dependency on `crypto-primitive`, we use their crate instead of
/// a copy-paste. We needed the newer `crypto-primitive` for serializing.
#[derive(Clone)]
pub(crate) struct IOPTranscript<F: PrimeField> {
    transcript: Transcript,
    is_empty: bool,
    #[doc(hidden)]
    phantom: PhantomData<F>,
}

// TODO: merge this with jf_plonk::transcript
impl<F: PrimeField> IOPTranscript<F> {
    /// Create a new IOP transcript.
    pub(crate) fn new(label: &'static [u8]) -> Self {
        Self {
            transcript: Transcript::new(label),
            is_empty: true,
            phantom: PhantomData,
        }
    }

    /// Append the message to the transcript.
    pub(crate) fn append_message(&mut self, label: &'static [u8], msg: &[u8]) -> Result<(), Error> {
        self.transcript.append_message(label, msg);
        self.is_empty = false;
        Ok(())
    }

    /// Append the message to the transcript.
    pub(crate) fn append_serializable_element<S: CanonicalSerialize>(
        &mut self,
        label: &'static [u8],
        group_elem: &S,
    ) -> Result<(), Error> {
        self.append_message(
            label,
            &to_bytes!(group_elem).map_err(|_| Error::TranscriptError)?,
        )
    }

    /// Generate the challenge from the current transcript
    /// and append it to the transcript.
    ///
    /// The output field element is statistical uniform as long
    /// as the field has a size less than 2^384.
    pub(crate) fn get_and_append_challenge(&mut self, label: &'static [u8]) -> Result<F, Error> {
        //  we need to reject when transcript is empty
        if self.is_empty {
            return Err(Error::TranscriptError);
        }

        let mut buf = [0u8; 64];
        self.transcript.challenge_bytes(label, &mut buf);
        let challenge = F::from_le_bytes_mod_order(&buf);
        self.append_serializable_element(label, &challenge)?;
        Ok(challenge)
    }
}

#[inline]
#[cfg(test)]
pub(crate) fn to_field<F: Field>(v: Vec<u64>) -> Vec<F> {
    v.iter().map(|x| F::from(*x)).collect::<Vec<F>>()
}

// TODO: replace by https://github.com/arkworks-rs/crypto-primitives/issues/112.
#[cfg(test)]
use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;

#[cfg(test)]
pub(crate) fn test_sponge<F: PrimeField>() -> PoseidonSponge<F> {
    use ark_crypto_primitives::sponge::{poseidon::PoseidonConfig, CryptographicSponge};
    use ark_std::test_rng;

    let full_rounds = 8;
    let partial_rounds = 31;
    let alpha = 17;

    let mds = vec![
        vec![F::one(), F::zero(), F::one()],
        vec![F::one(), F::one(), F::zero()],
        vec![F::zero(), F::one(), F::one()],
    ];

    let mut v = Vec::new();
    let mut ark_rng = test_rng();

    for _ in 0..(full_rounds + partial_rounds) {
        let mut res = Vec::new();

        for _ in 0..3 {
            res.push(F::rand(&mut ark_rng));
        }
        v.push(res);
    }
    let config = PoseidonConfig::new(full_rounds, partial_rounds, alpha, mds, v, 2, 1);
    PoseidonSponge::new(&config)
}
