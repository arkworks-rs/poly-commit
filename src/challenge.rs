use ark_ff::Field;
use ark_sponge::FieldBasedCryptographicSponge;

pub enum ChallengeStrategy {
    Multivariate,
    Univariate,
}

/// State stored for univariate generator
pub struct UnivariateGeneratorState<F: Field> {
    gen: F,
    next: F,
}

impl<F: Field> UnivariateGeneratorState<F> {
    fn new(gen: F) -> Self {
        Self { gen, next: gen }
    }

    fn get_next(&mut self) -> F {
        let result = self.next;
        self.next *= self.gen;
        result
    }
}

#[derive(Copy, Clone)]
pub enum ChallengeGenerator<'a, F: Field, S: 'a + FieldBasedCryptographicSponge<F>> {
    Multivariate(&'a mut S),
    Univariate(UnivariateGeneratorState<F>),
}

impl<'a, F: Field, S: 'a + FieldBasedCryptographicSponge<F>> ChallengeGenerator<'a, F, S> {
    /// Returns a challenge generator with multivariate strategy. Each challenge is freshly squeezed
    /// from a sponge.
    pub fn new_multivariate(sponge: &'a mut S) -> Self {
        Self::Multivariate(sponge)
    }

    /// Returns a challenge generator with univariate strategy. Each challenge is a power of one
    /// squeezed element from sponge.
    pub fn new_univariate(sponge: &mut S) -> Self {
        let gen = sponge.squeeze_native_field_elements(1)[0];
        let univariate_state = UnivariateGeneratorState { gen, next: gen };
        Self::Univariate(univariate_state)
    }

    /// Returns the next challenge generated.
    pub fn next_challenge(&mut self) -> F {
        if let Self::Multivariate(s) = &mut self {
            s.squeeze_native_field_elements(1)[0]
        } else if let Self::Univariate(s) = &mut self {
            s.get_next()
        }
    }
}
