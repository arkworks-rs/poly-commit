use ark_ff::PrimeField;
use ark_sponge::CryptographicSponge;

/// Challenge Generator (todo doc)
/// Note that mutable reference cannot be cloned.
#[derive(Clone)]
pub enum ChallengeGenerator<F: PrimeField, S: CryptographicSponge> {
    /// Each challenge is freshly squeezed from a sponge.
    Multivariate(S),
    /// Each challenge is a power of one squeezed element from sponge.
    ///
    /// `Univariate(generator, next_element)`
    Univariate(F, F),
}

impl<F: PrimeField, S: CryptographicSponge> ChallengeGenerator<F, S> {
    /// Returns a challenge generator with multivariate strategy. Each challenge is freshly squeezed
    /// from a sponge.
    pub fn new_multivariate(sponge: S) -> Self {
        Self::Multivariate(sponge)
    }

    /// Returns a challenge generator with univariate strategy. Each challenge is a power of one
    /// squeezed element from sponge.
    pub fn new_univariate(sponge: &mut S) -> Self {
        let gen = sponge.squeeze_field_elements(1)[0];
        Self::Univariate(gen, gen)
    }

    /// Returns the next challenge generated.
    pub fn next_challenge(&mut self) -> F {
        match self {
            Self::Multivariate(s) => s.squeeze_field_elements(1)[0],
            Self::Univariate(gen, next) => {
                let result = next.clone();
                *next *= *gen;
                result
            }
        }
    }

    /// Returns the next challenge generated where next challenge has `size` bits. Only works for
    /// multivariate generator.
    ///
    /// ## Panics
    /// This function will panic if `self` is univariate.
    pub fn next_challenge_of_size(&mut self, size: ark_sponge::FieldElementSize) -> F {
        match self {
            Self::Multivariate(s) => s.squeeze_field_elements_with_sizes(&[size])[0],
            _ => {
                panic!("`next_challenge_of_size` only supports multivariate generator.")
            }
        }
    }

    /// Returns the sponge state if `self` is multivariate.
    ///
    /// ## Panics
    /// This function will panic is `self` is univariate.
    pub fn into_sponge(self) -> S {
        match self {
            Self::Multivariate(s) => s,
            _ => panic!("only multivariate generator can be converted to sponge."),
        }
    }
}
