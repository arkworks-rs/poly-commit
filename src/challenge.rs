use ark_ff::PrimeField;
use ark_sponge::{CryptographicSponge, FieldElementSize};

/// `ChallengeGenerator` generates opening challenges using multivariate or univariate strategy.
/// For multivariate strategy, each challenge is freshly squeezed from a sponge.
/// For univariate strategy, each challenge is a power of one squeezed element from sponge.
///
/// Note that mutable reference cannot be cloned.
#[derive(Clone)]
pub enum ChallengeGenerator<F: PrimeField, S: CryptographicSponge> {
    /// Each challenge is freshly squeezed from a sponge.
    Multivariate(S, FieldElementSize),
    /// Each challenge is a power of one squeezed element from sponge.
    ///
    /// `Univariate(generator, next_element)`
    Univariate(F, F),
}

impl<F: PrimeField, S: CryptographicSponge> ChallengeGenerator<F, S> {
    /// Returns a challenge generator with multivariate strategy. Each challenge is freshly squeezed
    /// from a sponge.
    pub fn new_multivariate(sponge: S) -> Self {
        Self::new_multivariate_of_size(sponge, FieldElementSize::Full)
    }

    /// Returns a challenge generator with multivariate strategy.Each challenge is freshly squeezed
    /// from a sponge and has `size` bits.
    pub fn new_multivariate_of_size(sponge: S, size: FieldElementSize) -> Self {
        Self::Multivariate(sponge, size)
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
            // multivariate (full)
            Self::Multivariate(sponge, FieldElementSize::Full) => {
                sponge.squeeze_field_elements(1)[0]
            }
            // multivariate (truncated)
            Self::Multivariate(sponge, size) => {
                sponge.squeeze_field_elements_with_sizes(&[*size])[0]
            }
            // univariate
            Self::Univariate(gen, next) => {
                let result = next.clone();
                *next *= *gen;
                result
            }
        }
    }

    /// Returns the sponge state if `self` is multivariate.
    ///
    /// ## Panics
    /// This function will panic if `self` is univariate.
    pub fn into_sponge(self) -> S {
        match self {
            Self::Multivariate(s, _) => s,
            _ => panic!("only multivariate generator can be converted to sponge."),
        }
    }
}
