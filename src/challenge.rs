use ark_ff::PrimeField;
use ark_sponge::FieldBasedCryptographicSponge;

/// State stored for univariate generator
pub struct UnivariateGeneratorState<F: PrimeField> {
    gen: F,
    next: F,
}

impl<F: PrimeField> UnivariateGeneratorState<F> {
    fn new(gen: F) -> Self {
        Self { gen, next: gen }
    }

    fn get_next(&mut self) -> F {
        let result = self.next;
        self.next *= self.gen;
        result
    }
}

// TODO: Copy and Clone trait cannot be derived
/// Challenge Generator (todo doc)
pub enum ChallengeGenerator<'a, F: PrimeField, S: 'a + FieldBasedCryptographicSponge<F>> {
    /// todo: doc
    Multivariate(&'a mut S),
    /// todo: doc
    Univariate(UnivariateGeneratorState<F>),
}

impl<'a, F: PrimeField, S: 'a + FieldBasedCryptographicSponge<F>> ChallengeGenerator<'a, F, S> {
    /// Returns a challenge generator with multivariate strategy. Each challenge is freshly squeezed
    /// from a sponge.
    pub fn new_multivariate(sponge: &'a mut S) -> Self {
        Self::Multivariate(sponge)
    }

    /// Returns a challenge generator with univariate strategy. Each challenge is a power of one
    /// squeezed element from sponge.
    pub fn new_univariate(sponge: &mut S) -> Self {
        let gen = sponge.squeeze_native_field_elements(1)[0];
        let univariate_state = UnivariateGeneratorState::new(gen);
        Self::Univariate(univariate_state)
    }

    /// Returns the next challenge generated.
    pub fn next_challenge(&mut self) -> F {
        if let Self::Multivariate(s) = self {
            s.squeeze_native_field_elements(1)[0]
        } else if let Self::Univariate(s) = self {
            s.get_next()
        } else {
            // should not happen
            panic!()
        }
    }
}
