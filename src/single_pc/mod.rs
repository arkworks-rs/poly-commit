use crate::*;
use algebra::Field;
use rand_core::RngCore;

/// Describes the interface for a polynomial commitment scheme that allows
/// a sender to commit to a single polynomial and later provide a succinct proof
/// of evaluation for that commitment.
/// ```compile_fail
///
/// let rng = &mut thread_rng();
/// // Generate committer and verifier keys for a maximum polynomial degree
/// // of 10.
/// let degree = 10;
/// let (ck, vk) = SinglePC::setup(degree, rng);
///
/// // Generate a random degree-10 polynomial.
/// let p = Polynomial::rand(degree, rng);;
/// let hiding_bound = Some(1);
///
/// // Commit to this polynomial
/// let (comm, rand) = SinglePC::commit(&ck, &p, hiding_bound, rng)?;
///
/// // Evaluate the polynomial at a random point, and generate an evaluation
/// // proof.
/// let point = F::rand(rng);
/// let value = p.evaluate(point);
/// let proof = SinglePC::open(&ck, &p, point, &rand)?;
///
/// // Verify the evaluation proof.
/// assert!(SinglePC::check(&vk, &comm, point, value, &proof)?, "proof was incorrect");
/// ```
pub trait SinglePolynomialCommitment<F: Field> {
    /// The universal parameters that must be specialized for specific degrees.
    type UniversalParams: PCUniversalParams;
    /// The committer key for the scheme; used to commit to a polynomial and then
    /// open the commitment to produce an evaluation proof.
    type CommitterKey: PCCommitterKey;
    /// The verifier key for the scheme; used to check an evaluation proof.
    type VerifierKey: PCVerifierKey;
    /// The commitment to a polynomial.
    type Commitment: PCCommitment;
    /// The commitment randomness.
    type Randomness: PCRandomness;
    /// The evaluation proof.
    type Proof: Clone;
    /// The evaluation proof for a query set.
    type BatchProof: Clone + From<Vec<Self::Proof>>;
    /// The error type for the scheme.
    type Error: std::error::Error;

    /// Constructs universal public parameters when given as input the maximum 
    /// degree `max_degree` for the polynomial commitment scheme.
    fn setup<R: RngCore>(
        max_degree: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error>;

    /// "Trims" the universal parameters to support polynomials having coefficients
    /// defined by `coefficient_support`.
    fn trim(
        pp: &Self::UniversalParams,
        coefficient_support: CoefficientSupport,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error>;

    /// Outputs a commitment to `polynomial`. If `hiding_bound.is_some()`, then the
    /// resulting commitment is hiding up to `hiding_bound.unwrap()` number of queries.
    /// `rng` should not be `None` if `hiding_bound.is_some()`.
    /// If `hiding_bound.is_none()`, then the randomness is `Self::Randomness::empty()`.
    fn commit(
        ck: &Self::CommitterKey,
        polynomial: &Polynomial<F>,
        hiding_bound: Option<usize>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<(Self::Commitment, Self::Randomness), Self::Error>;

    /// On input a polynomial `p` and a point `point` in the field `F`,
    /// `open` outputs an evaluation proof for `p` at `point`.
    fn open(
        ck: &Self::CommitterKey,
        p: &Polynomial<F>,
        point: F,
        r: &Self::Randomness,
    ) -> Result<Self::Proof, Self::Error>;

    /// Verifies that `value` is the evaluation at `point` of the polynomial
    /// committed inside `comm`.
    fn check(
        vk: &Self::VerifierKey,
        comm: &Self::Commitment,
        point: F,
        value: F,
        proof: &Self::Proof,
    ) -> Result<bool, Self::Error>;



    /// Check a batch of proofs
    fn batch_check<R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: &[Self::Commitment],
        points: &[F],
        values: &[F],
        proofs: &[Self::Proof],
        rng: &mut R,
    ) -> Result<bool, Self::Error>;
}

/// Implements the [KZG10](kzg10) construction that satisfies the `SinglePolynomialCommitment`
/// trait.
///
/// [kzg10]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
pub mod kzg10;

/// A struct that represents the non-zero coefficients of a polynomial.
pub(crate) struct CoefficientSupport {
    lower_interval_bound: usize,
    upper_interval_bound: usize,
}

impl CoefficientSupport {
    fn from_dual_interval(lower_interval_bound: usize, upper_interval_bound: usize) -> Self {
        Self {
            lower_interval_bound,
            upper_interval_bound,
        }
    }
}

#[cfg(test)]
pub mod tests {
    use crate::*;
    use algebra::Field;
    use algebra::UniformRand;
    use rand::thread_rng;

    pub fn end_to_end_test<F, SinglePC>() -> Result<(), SinglePC::Error>
    where
        F: Field,
        SinglePC: SinglePolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        for _ in 0..100 {
            let mut degree = 0;
            while degree <= 1 {
                degree = usize::rand(rng) % 20;
            }
            let pp = SinglePC::setup(degree, rng)?;
            let bounds = 
            let (ck, vk) = SinglePC::trim(degree, degree);
            let p = Polynomial::rand(degree, rng);;
            let hiding_bound = Some(1);
            let (comm, rand) = SinglePC::commit(&ck, &p, hiding_bound, Some(rng))?;
            let point = F::rand(rng);
            let value = p.evaluate(point);
            let proof = SinglePC::open(&ck, &p, point, &rand)?;
            assert!(
                SinglePC::check(&vk, &comm, point, value, &proof)?,
                "proof was incorrect for max_degree = {}, polynomial_degree = {}, hiding_bound = {:?}",
                degree,
                p.degree(),
                hiding_bound,
            );
        }
        Ok(())
    }

    pub fn linear_polynomial_test<F, SinglePC>() -> Result<(), SinglePC::Error>
    where
        F: Field,
        SinglePC: SinglePolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        for _ in 0..100 {
            let degree = 2;
            let (ck, vk) = SinglePC::setup(degree, rng)?;
            let p = Polynomial::rand(1, rng);;
            let hiding_bound = Some(1);
            let (comm, rand) = SinglePC::commit(&ck, &p, hiding_bound, Some(rng))?;
            let point = F::rand(rng);
            let value = p.evaluate(point);
            let proof = SinglePC::open(&ck, &p, point, &rand)?;
            assert!(
                SinglePC::check(&vk, &comm, point, value, &proof)?,
                "proof was incorrect for max_degree = {}, polynomial_degree = {}, hiding_bound = {:?}",
                degree,
                p.degree(),
                hiding_bound,
            );
        }
        Ok(())
    }

    pub fn batch_check_test<F, SinglePC>() -> Result<(), SinglePC::Error>
    where
        F: Field,
        SinglePC: SinglePolynomialCommitment<F>,
    {
        let rng = &mut thread_rng();
        for _ in 0..10 {
            let mut degree = 0;
            while degree <= 1 {
                degree = usize::rand(rng) % 20;
            }
            let (ck, vk) = SinglePC::setup(degree, rng)?;
            let mut comms = Vec::new();
            let mut values = Vec::new();
            let mut points = Vec::new();
            let mut proofs = Vec::new();
            for _ in 0..10 {
                let p = Polynomial::rand(degree, rng);;
                let hiding_bound = Some(1);
                let (comm, rand) = SinglePC::commit(&ck, &p, hiding_bound, Some(rng))?;
                let point = F::rand(rng);
                let value = p.evaluate(point);
                let proof = SinglePC::open(&ck, &p, point, &rand)?;

                assert!(SinglePC::check(&vk, &comm, point, value, &proof)?);
                comms.push(comm);
                values.push(value);
                points.push(point);
                proofs.push(proof);
            }
            assert!(SinglePC::batch_check(
                &vk, &comms, &points, &values, &proofs, rng
            )?);
        }
        Ok(())
    }
}
