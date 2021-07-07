use ark_ff::PrimeField;
use ark_sponge::CryptographicSponge;

/// Challenge Generator (todo doc)
/// TODO: probably move it to sponge
/// Note that mutable reference cannot be cloned.
pub enum ChallengeGenerator<'a, F: PrimeField, S: 'a + CryptographicSponge> {
    /// Each challenge is freshly squeezed from a sponge.
    Multivariate(&'a mut S),
    /// Each challenge is a power of one squeezed element from sponge.
    ///
    /// `Univariate(generator, next_element)`
    Univariate(F, F),
}

impl<'a, F: PrimeField, S: 'a + CryptographicSponge> ChallengeGenerator<'a, F, S> {
    /// Returns a challenge generator with multivariate strategy. Each challenge is freshly squeezed
    /// from a sponge.
    pub fn new_multivariate(sponge: &'a mut S) -> Self {
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

    // TODO: pub fn next_challenge_with_bit_size -> Option<F>
}
