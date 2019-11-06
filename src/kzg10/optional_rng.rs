use rand_core::RngCore;

// This trick is necessary because `Option<&mut R>` is not implicitly reborrowed
// like `&mut R` is. As a result, we define a dummy rng here that should be used
// when `commit` gets `rng = None`
//
// Basically, we define a "dummy rng" that does nothing
// (corresponding to the case that `rng = None`).
pub(crate) struct OptionalRng<R>(pub(crate) Option<R>);

impl<R: RngCore> RngCore for OptionalRng<R> {
    #[inline]
    fn next_u32(&mut self) -> u32 {
        (&mut self.0).as_mut().map_or(0, |r| r.next_u32())
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        (&mut self.0).as_mut().map_or(0, |r| r.next_u64())
    }

    #[inline]
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        (&mut self.0).as_mut().map_or((), |r| r.fill_bytes(dest))
    }

    #[inline]
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        Ok(self.fill_bytes(dest))
    }
}
