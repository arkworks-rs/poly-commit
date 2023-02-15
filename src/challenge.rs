use ark_crypto_primitives::sponge::{CryptographicSponge, FieldElementSize};
use ark_ff::PrimeField;

/// `ChallengeGenerator` generates opening challenges using multivariate or univariate strategy.
/// For multivariate strategy, each challenge is freshly squeezed from a sponge.
/// For univariate strategy, each challenge is a power of one squeezed element from sponge.
///
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

    /// Returns a challenge of size `size`.
    /// * If `self == Self::Multivariate(...)`, then this squeezes out a challenge of size `size`.
    /// * If `self == Self::Univariate(...)`, then this ignores the `size` argument and simply squeezes out
    /// the next field element.
    pub fn try_next_challenge_of_size(&mut self, size: FieldElementSize) -> F {
        match self {
            // multivariate (full)
            Self::Multivariate(sponge) => sponge.squeeze_field_elements_with_sizes(&[size])[0],
            // univariate
            Self::Univariate(gen, next) => {
                let result = next.clone();
                *next *= *gen;
                result
            }
        }
    }
    /// Returns the next challenge generated.
    pub fn next_challenge(&mut self) -> F {
        self.try_next_challenge_of_size(FieldElementSize::Full)
    }

    /// Returns the sponge state if `self` is multivariate. Returns `None` otherwise.
    pub fn into_sponge(self) -> Option<S> {
        match self {
            Self::Multivariate(s) => Some(s),
            _ => None,
        }
    }
}
