use ark_crypto_primitives::sponge::{CryptographicSponge, FieldElementSize};
use ark_ff::PrimeField;

/// `ChallengeGenerator` generates opening challenges using independent or correlated strategy.
/// For independent strategy, each challenge is freshly squeezed from a sponge.
/// For correlated strategy, each challenge is a power of one squeezed element from sponge.
///
/// Note that mutable reference cannot be cloned.
#[derive(Clone)]
pub enum ChallengeGenerator<F: PrimeField, S: CryptographicSponge> {
    /// Each challenge is freshly squeezed from a sponge.
    Independent(S),
    /// Each challenge is a power of one squeezed element from sponge.
    ///
    /// `Correlated(generator, next_element)`
    Correlated(F, F),
}

impl<F: PrimeField, S: CryptographicSponge> ChallengeGenerator<F, S> {
    /// Returns a challenge generator with independent strategy. Each challenge is freshly squeezed
    /// from a sponge.
    pub fn new_independent(sponge: S) -> Self {
        Self::Independent(sponge)
    }

    /// Returns a challenge generator with correlated strategy. Each challenge is a power of one
    /// squeezed element from sponge.
    pub fn new_correlated(sponge: &mut S) -> Self {
        let gen = sponge.squeeze_field_elements(1)[0];
        Self::Correlated(gen, gen)
    }

    /// Returns a challenge of size `size`.
    /// * If `self == Self::Independent(...)`, then this squeezes out a challenge of size `size`.
    /// * If `self == Self::Correlated(...)`, then this ignores the `size` argument and simply squeezes out
    /// the next field element.
    pub fn try_next_challenge_of_size(&mut self, size: FieldElementSize) -> F {
        match self {
            // independent (full)
            Self::Independent(sponge) => sponge.squeeze_field_elements_with_sizes(&[size])[0],
            // correlated
            Self::Correlated(gen, next) => {
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

    /// Returns the sponge state if `self` is independent. Returns `None` otherwise.
    pub fn into_sponge(self) -> Option<S> {
        match self {
            Self::Independent(s) => Some(s),
            _ => None,
        }
    }
}
