use crate::multi_pc::{Evaluations, QuerySet};
use crate::*;
use derivative::Derivative;
use std::collections::{BTreeSet, BTreeMap};
use std::fmt::Debug;
use std::hash::Hash;
use std::marker::PhantomData;

use algebra::PrimeField;
use rand::RngCore;

/// Generic construction of a `MultiPolynomialCommitment` scheme from a
/// `SinglePolynomialCommitment` scheme whenever the commitment and randomness of the
/// `SinglePolynomialCommitment` scheme are additively homomorphic.
/// Specifically, we require `C = MultiPolynomialCommitment::Commitment`
///
/// The construction follows the blueprint laid out in [CHMMVW19](insert eprint link).
pub struct MultiPCFromSinglePC<F, SinglePC> {
    _field: PhantomData<F>,
    _spc: PhantomData<SinglePC>,
}

/// A trait that signals that `Self::Commitment` can be combined efficiently.
pub trait SinglePCExt<F: Field>: SinglePolynomialCommitment<F> {
    /// Take a linear combination of `commitments`.
    fn combine_commitments(commitments: &[Self::Commitment], coeffs: &[F]) -> Self::Commitment;
}

/// Commitment to a polynomial that optionally enforces a degree bound.
/// Output by `MultiPCFromSinglePC::commit`.
#[derive(Derivative)]
#[derivative(
    Default(bound = "SinglePC::Commitment: Default"),
    Hash(bound = "SinglePC::Commitment: Hash"),
    Clone(bound = "SinglePC::Commitment: Clone"),
    Copy(bound = "SinglePC::Commitment: Copy"),
    Debug(bound = "SinglePC::Commitment: Debug"),
    PartialEq(bound = "SinglePC::Commitment: PartialEq"),
    Eq(bound = "SinglePC::Commitment: Eq")
)]
pub struct Commitment<F: PrimeField, SinglePC: SinglePolynomialCommitment<F>> {
    comm: SinglePC::Commitment,
    shifted_comm: Option<SinglePC::Commitment>,
}

impl<F: PrimeField, SinglePC: SinglePolynomialCommitment<F>> PCCommitment
    for Commitment<F, SinglePC>
{
    fn empty() -> Self {
        Self {
            comm: SinglePC::Commitment::empty(),
            shifted_comm: Some(SinglePC::Commitment::empty()),
        }
    }

    fn has_degree_bound(&self) -> bool {
        self.shifted_comm.is_some()
    }

    fn size_in_bytes(&self) -> usize {
        self.comm.size_in_bytes() + self.shifted_comm.as_ref().map_or(0, |c| c.size_in_bytes())
    }
}

impl<F: PrimeField, SinglePC: SinglePolynomialCommitment<F>> algebra::ToBytes
    for Commitment<F, SinglePC>
{
    #[inline]
    fn write<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        self.comm.write(&mut writer)?;
        let shifted_exists = self.shifted_comm.is_some();
        shifted_exists.write(&mut writer)?;
        self.shifted_comm
            .as_ref()
            .unwrap_or(&SinglePC::Commitment::empty())
            .write(&mut writer)
    }
}

/// Randomness used to make `Commitment` hiding. Output by `MultiPCFromSinglePC::commit`.
#[derive(Derivative)]
#[derivative(
    Default(bound = "SinglePC::Randomness: Default"),
    Hash(bound = "SinglePC::Randomness: Hash"),
    Clone(bound = "SinglePC::Randomness: Clone"),
    Copy(bound = "SinglePC::Randomness: Copy"),
    Debug(bound = "SinglePC::Randomness: Debug"),
    PartialEq(bound = "SinglePC::Randomness: PartialEq"),
    Eq(bound = "SinglePC::Randomness: Eq")
)]
pub struct Randomness<F: PrimeField, SinglePC: SinglePolynomialCommitment<F>> {
    rand: SinglePC::Randomness,
    shifted_rand: Option<SinglePC::Randomness>,
}

impl<F: PrimeField, SinglePC: SinglePolynomialCommitment<F>> PCRandomness
    for Randomness<F, SinglePC>
{
    fn empty() -> Self {
        Self {
            rand: SinglePC::Randomness::empty(),
            shifted_rand: Some(SinglePC::Randomness::empty()),
        }
    }

    fn rand<R: RngCore>(hiding_bound: usize, rng: &mut R) -> Self {
        Self {
            rand: SinglePC::Randomness::rand(hiding_bound, rng),
            shifted_rand: Some(SinglePC::Randomness::rand(hiding_bound, rng)),
        }
    }
}

/// Evaluation proof output by `MultiPCFromSinglePC::open`.
#[derive(Derivative)]
#[derivative(
    Default(bound = "SinglePC::Proof: Default"),
    Hash(bound = "SinglePC::Proof: Hash"),
    Clone(bound = "SinglePC::Proof: Clone"),
    Debug(bound = "SinglePC::Proof: Debug"),
    PartialEq(bound = "SinglePC::Proof: PartialEq"),
    Eq(bound = "SinglePC::Proof: Eq")
)]
pub struct Proof<F: PrimeField, SinglePC: SinglePolynomialCommitment<F>> {
    proofs: Vec<SinglePC::Proof>,
}

pub(crate) fn shift_polynomial<F: Field>(
    p: &Polynomial<F>,
    degree_bound: usize,
    max_degree: usize,
) -> Polynomial<F> {
    if p.is_zero() {
        Polynomial::zero()
    } else {
        let mut shifted_polynomial_coeffs = vec![F::zero(); max_degree - degree_bound];
        shifted_polynomial_coeffs.extend_from_slice(&p.coeffs);
        Polynomial::from_coefficients_vec(shifted_polynomial_coeffs)
    }
}

/// Error type for `MultiPCFromSinglePC`.
#[derive(Debug)]
pub enum Error<E> {
    /// The degree of the `index`-th polynomial passed to `commit` or `open`
    /// was too large.
    PolynomialDegreeTooLarge {
        /// Degree of the polynomial.
        poly_degree: usize,
        /// Maximum supported degree.
        max_degree: usize,
        /// Index of the offending polynomial.
        label: String,
    },
    /// The degree bound for the `index`-th polynomial passed to `commit`, `open`
    /// or `check` was incorrect, that is, `degree_bound >= poly_degree` or
    /// `degree_bound <= max_degree`.
    IncorrectDegreeBound {
        /// Degree of the polynomial.
        poly_degree: usize,
        /// Degree bound.
        degree_bound: usize,
        /// Maximum supported degree.
        max_degree: usize,
        /// Index of the offending polynomial.
        label: String,
    },
    /// The inputs to `commit`, `open` or `verify` had incorrect lengths.
    IncorrectInputLength(String),
    /// The query set referenced a non-existent polynomial.
    IncorrectQuerySet(&'static str),
    /// The evaluations referenced a non-existent member of the query set.
    IncorrectEvaluation(&'static str),
    /// An error from the underlying `SinglePC`.
    SPCError(E),
}

impl<E> From<E> for Error<E> {
    fn from(other: E) -> Self {
        Error::SPCError(other)
    }
}

impl<E> Error<E> {
    fn poly_degree_too_large(poly_degree: usize, max_degree: usize, label: String) -> Self {
        Error::PolynomialDegreeTooLarge {
            poly_degree,
            max_degree,
            label,
        }
    }

    fn incorrect_bound(
        poly_degree: usize,
        degree_bound: usize,
        max_degree: usize,
        label: String,
    ) -> Self {
        Error::IncorrectDegreeBound {
            poly_degree,
            degree_bound,
            max_degree,
            label,
        }
    }

    fn check_degrees(
        d: usize,
        bound: Option<usize>,
        max_degree: usize,
        label: String,
    ) -> Result<(), Self> {
        if let Some(bound) = bound {
            if d > max_degree {
                Err(Error::poly_degree_too_large(d, max_degree, label))
            } else if bound < d || bound > max_degree {
                Err(Error::incorrect_bound(d, bound, max_degree, label))
            } else {
                Ok(())
            }
        } else {
            Ok(())
        }
    }
}

impl<E: std::fmt::Display> std::fmt::Display for Error<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::PolynomialDegreeTooLarge {
                poly_degree,
                max_degree,
                label,
            } => write!(
                f,
                "the degree of the polynomial {} ({:?}) is greater than\
                 the maximum supported degree ({:?})",
                label, poly_degree, max_degree
            ),
            Error::IncorrectDegreeBound {
                poly_degree,
                degree_bound,
                max_degree,
                label,
            } => write!(
                f,
                "the degree bound ({:?}) for the polynomial {} \
                 (having degree {:?}) is greater than the maximum \
                 supported degree ({:?})",
                degree_bound, label, poly_degree, max_degree
            ),
            Error::IncorrectInputLength(err) => write!(f, "{}", err),
            Error::IncorrectQuerySet(err) => write!(f, "{}", err),
            Error::IncorrectEvaluation(err) => write!(f, "{}", err),
            Error::SPCError(err) => write!(f, "{}", err),
        }
    }
}

impl<E: std::error::Error> std::error::Error for Error<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

impl<F, SinglePC> MultiPolynomialCommitment<F> for MultiPCFromSinglePC<F, SinglePC>
where
    F: PrimeField,
    SinglePC: SinglePolynomialCommitment<F> + SinglePCExt<F>,
    for<'a> SinglePC::Commitment: std::ops::AddAssign<(F, &'a SinglePC::Commitment)>,
    for<'a> SinglePC::Randomness: std::ops::AddAssign<(F, &'a SinglePC::Randomness)>,
{
    type CommitterKey = SinglePC::CommitterKey;
    type VerifierKey = SinglePC::VerifierKey;
    type Commitment = Commitment<F, SinglePC>;
    type Randomness = Randomness<F, SinglePC>;
    type Proof = Proof<F, SinglePC>;
    type Error = Error<SinglePC::Error>;

    /// Constructs public parameters when given as input the maximum degree `degree`
    /// for the polynomial commitment scheme.
    fn setup<R: RngCore>(
        degree: usize,
        rng: &mut R,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        Ok(SinglePC::setup(degree, rng)?)
    }

    /// Outputs a commitment to `polynomial`.
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, F>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<(Vec<LabeledCommitment<Self::Commitment>>, Vec<Self::Randomness>), Self::Error> {
        let commit_time = start_timer!(|| "Committing to polynomials");

        let mut commitments = Vec::new();
        let mut randomness = Vec::new();
        let max_degree = ck.max_degree();

        let rng = &mut optional_rng::OptionalRng(rng);
        for polynomial in polynomials {
            let label = polynomial.label();
            let degree_bound = polynomial.degree_bound();
            let hiding_bound = polynomial.hiding_bound();
            let polynomial = polynomial.polynomial();

            Error::check_degrees(polynomial.degree(), degree_bound, max_degree, label.to_string())?;

            let commit_time = start_timer!(|| format!(
                "{} of degree {}, bound {:?}, and hiding bound {:?}",
                _label,
                polynomial.degree(),
                degree_bound,
                hiding_bound,
            ));
            let (comm, rand) = SinglePC::commit(ck, polynomial, hiding_bound, Some(rng))?;
            let (shifted_comm, shifted_rand) = if let Some(degree_bound) = degree_bound {
                if degree_bound < max_degree {
                    let s_polynomial = shift_polynomial(polynomial, degree_bound, max_degree);
                    assert!(
                        polynomial.degree() <= s_polynomial.degree()
                            && s_polynomial.degree() <= max_degree
                            && s_polynomial.degree()
                                == polynomial.degree() + max_degree - degree_bound,
                        "polynomial.degree(): {}; s_polynomial.degree(): {}; max_degree: {}.",
                        polynomial.degree(),
                        s_polynomial.degree(),
                        max_degree
                    );
                    let (shifted_comm, shifted_rand) =
                        SinglePC::commit(ck, &s_polynomial, hiding_bound, Some(rng))?;
                    (Some(shifted_comm), Some(shifted_rand))
                } else {
                    (None, None)
                }
            } else {
                (None, None)
            };

            let comm = Commitment { comm, shifted_comm };
            let rand = Randomness { rand, shifted_rand };
            commitments.push(LabeledCommitment::new(label.clone(), comm, degree_bound));
            randomness.push(rand);
            end_timer!(commit_time);
        }
        end_timer!(commit_time);
        Ok((commitments, randomness))
    }

    /// On input a polynomial `p` and a point `point`, outputs `p(point)` and a proof for the same.
    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, F>>,
        query_set: &QuerySet<F>,
        opening_challenge: F,
        rands: &[Self::Randomness],
    ) -> Result<Self::Proof, Self::Error> {
        let open_time = start_timer!(|| format!(
            "Opening {} polynomials at query set of size {}",
            polynomials.len(),
            query_set.len(),
        ));

        let (polynomials, rands): (BTreeMap<_, _>, BTreeMap<_, _>) = labeled_polynomials
            .into_iter()
            .zip(rands)
            .map(|(poly, r)| ((poly.label(), poly), (poly.label(), r)))
            .unzip();

        let mut query_to_labels_map = BTreeMap::new();

        for (label, point) in query_set.iter() {
            let labels = query_to_labels_map.entry(point).or_insert(BTreeSet::new());
            labels.insert(label);
        }

        let mut proofs = Vec::new();
        let max_degree = ck.max_degree();
        for (query, labels) in query_to_labels_map.into_iter() {
            let mut p = Polynomial::zero();
            let mut r = SinglePC::Randomness::empty();
            let lc_time =
                start_timer!(|| format!("Randomly combining {} polynomials", labels.len()));
            for (j, label) in labels.into_iter().enumerate() {
                let polynomial = polynomials.get(label).ok_or(Error::IncorrectQuerySet(
                        "query set references polynomial with incorrect label",
                    ))?;

                let rand = rands.get(label).ok_or(Error::IncorrectQuerySet(
                        "query set references randomness with incorrect label",
                    ))?;
                let degree_bound = polynomial.degree_bound();

                Error::check_degrees(polynomial.degree(), degree_bound, max_degree, label.to_string())?;

                // compute challenge^j and challenge^{j+1}.
                let challenge_j = opening_challenge.pow([2 * j as u64]);

                assert_eq!(
                    degree_bound.is_some(),
                    rand.shifted_rand.is_some()
                );

                p += (challenge_j, polynomial.polynomial());
                r += (challenge_j, &rand.rand);

                if let Some(degree_bound) = degree_bound {
                    let challenge_j_1 = challenge_j * &opening_challenge;

                    let s_polynomial = shift_polynomial(&polynomial, degree_bound, max_degree);

                    p += (challenge_j_1, &s_polynomial);
                    r += (challenge_j_1, &rand.shifted_rand.as_ref().unwrap());
                }
            }
            end_timer!(lc_time);
            let proof_time = start_timer!(|| "Creating SinglePC::Proof");
            let proof = SinglePC::open(ck, &p, *query, &r)?;
            end_timer!(proof_time);

            proofs.push(proof);
        }
        end_timer!(open_time);

        Ok(Proof { proofs })
    }

    /// Verifies that `value` is truly the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    fn check<R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: &[LabeledCommitment<Self::Commitment>],
        query_set: &QuerySet<F>,
        values: &Evaluations<F>,
        proof: &Self::Proof,
        opening_challenge: F,
        rng: &mut R,
    ) -> Result<bool, Self::Error> {
        let commitments: BTreeMap<_, _> = commitments.into_iter().map(|c| (c.label(), c)).collect();
        let mut query_to_labels_map = BTreeMap::new();

        for (label, point) in query_set.iter() {
            let labels = query_to_labels_map.entry(point).or_insert(BTreeSet::new());
            labels.insert(label);
        }
        assert_eq!(proof.proofs.len(), query_to_labels_map.len());

        let max_degree = vk.max_degree();
        let mut combined_comms = Vec::new();
        let mut combined_queries = Vec::new();
        let mut combined_evals = Vec::new();
        for (query, labels) in query_to_labels_map.into_iter() {
            let lc_time =
                start_timer!(|| format!("Randomly combining {} commitments", indices.len()));
            let mut comms_to_combine = Vec::new();
            let mut values_to_combine = Vec::new();
            let mut randomizers = Vec::new();
            let mut challenge_j = F::one();
            for label in labels.into_iter() {
                let commitment = commitments.get(label).ok_or(Error::IncorrectQuerySet(
                        "query set references commitment with incorrect label",
                    ))?;
                let degree_bound = commitment.degree_bound();
                let commitment = commitment.commitment();
                assert_eq!(
                    degree_bound.is_some(),
                    commitment.shifted_comm.is_some()
                );

                let v_i = values
                    .get(&(label.clone(), *query))
                    .ok_or(Error::IncorrectEvaluation("value does not exist"))?;

                comms_to_combine.push(commitment.comm.clone());
                values_to_combine.push(*v_i);
                randomizers.push(challenge_j);

                if let Some(degree_bound) = degree_bound {
                    let challenge_j_1 = challenge_j * &opening_challenge;
                    let shift = query.pow([(max_degree - degree_bound) as u64]);

                    comms_to_combine.push(commitment.shifted_comm.as_ref().unwrap().clone());
                    values_to_combine.push(shift * v_i);
                    randomizers.push(challenge_j_1);
                }
                challenge_j *= &opening_challenge.square();
            }
            let v = values_to_combine
                .into_iter()
                .zip(&randomizers)
                .fold(F::zero(), |x, (v, r)| x + &(v * r));
            let c = SinglePC::combine_commitments(&comms_to_combine, &randomizers);

            end_timer!(lc_time);
            combined_comms.push(c);
            combined_queries.push(*query);
            combined_evals.push(v);
        }
        let proof_time = start_timer!(|| "Checking SinglePC::Proof");
        let result = SinglePC::batch_check(
            vk,
            &combined_comms,
            &combined_queries,
            &combined_evals,
            &proof.proofs,
            rng,
        )?;
        end_timer!(proof_time);
        Ok(result)
    }
}

// This trick is necessary because `Option<&mut R>` is not implicitly reborrowed
// like `&mut R` is. As a result, we define a dummy rng here that should be used
// when `commit` gets `rng = None`
//
// Basically, we define a "dummy rng" that does nothing
// (corresponding to the case that `rng = None`).
pub(super) mod optional_rng {
    use rand::RngCore;
    pub(super) struct OptionalRng<R>(pub(super) Option<R>);

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
        fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
            Ok(self.fill_bytes(dest))
        }
    }
}

mod impl_kzg10 {
    use super::*;
    use crate::single_pc::kzg10::*;
    use algebra::{AffineCurve, PairingEngine, ProjectiveCurve};
    impl<E: PairingEngine> SinglePCExt<E::Fr> for KZG10<E> {
        fn combine_commitments(comms: &[Self::Commitment], coeffs: &[E::Fr]) -> Self::Commitment {
            let mut result = E::G1Projective::zero();
            for (comm, coeff) in comms.iter().zip(coeffs) {
                result += &comm.0.mul(*coeff);
            }
            Commitment(result.into())
        }
    }
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]

    use crate::multi_pc::mpc_from_spc::*;
    use crate::single_pc::kzg10::KZG10;
    use algebra::curves::bls12_377::Bls12_377;
    use algebra::curves::bls12_381::Bls12_381;
    use algebra::curves::mnt6::MNT6;
    use algebra::curves::sw6::SW6;
    use algebra::PairingEngine;

    type MultiPC<E> = MultiPCFromSinglePC<<E as PairingEngine>::Fr, KZG10<E>>;
    type MultiPC_Bls12_381 = MultiPC<Bls12_381>;
    type MultiPC_Bls12_377 = MultiPC<Bls12_377>;
    type MultiPC_MNT6 = MultiPC<MNT6>;
    type MultiPC_SW6 = MultiPC<SW6>;

    #[test]
    fn single_poly_test() {
        use crate::multi_pc::tests::*;
        single_poly_test::<_, MultiPC_Bls12_377>().expect("test failed for bls12-377");
        single_poly_test::<_, MultiPC_Bls12_381>().expect("test failed for bls12-381");
        single_poly_test::<_, MultiPC_MNT6>().expect("test failed for MNT6");
        single_poly_test::<_, MultiPC_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn single_poly_degree_bound_test() {
        use crate::multi_pc::tests::*;
        single_poly_degree_bound_test::<_, MultiPC_Bls12_377>().expect("test failed for bls12-377");
        single_poly_degree_bound_test::<_, MultiPC_Bls12_381>().expect("test failed for bls12-381");
        single_poly_degree_bound_test::<_, MultiPC_MNT6>().expect("test failed for MNT6");
        single_poly_degree_bound_test::<_, MultiPC_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn single_poly_degree_bound_multiple_queries_test() {
        use crate::multi_pc::tests::*;
        single_poly_degree_bound_multiple_queries_test::<_, MultiPC_Bls12_377>()
            .expect("test failed for bls12-377");
        single_poly_degree_bound_multiple_queries_test::<_, MultiPC_Bls12_381>()
            .expect("test failed for bls12-381");
        single_poly_degree_bound_multiple_queries_test::<_, MultiPC_MNT6>()
            .expect("test failed for MNT6");
        single_poly_degree_bound_multiple_queries_test::<_, MultiPC_SW6>()
            .expect("test failed for SW6");
    }

    #[test]
    fn two_polys_degree_bound_single_query_test() {
        use crate::multi_pc::tests::*;
        two_polys_degree_bound_single_query_test::<_, MultiPC_Bls12_377>()
            .expect("test failed for bls12-377");
        two_polys_degree_bound_single_query_test::<_, MultiPC_Bls12_381>()
            .expect("test failed for bls12-381");
        two_polys_degree_bound_single_query_test::<_, MultiPC_MNT6>()
            .expect("test failed for MNT6");
        two_polys_degree_bound_single_query_test::<_, MultiPC_SW6>().expect("test failed for SW6");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::multi_pc::tests::*;
        full_end_to_end_test::<_, MultiPC_Bls12_377>().expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_test::<_, MultiPC_Bls12_381>().expect("test failed for bls12-381");
        println!("Finished bls12-381");
        full_end_to_end_test::<_, MultiPC_MNT6>().expect("test failed for MNT6");
        println!("Finished mnt6");
        full_end_to_end_test::<_, MultiPC_SW6>().expect("test failed for SW6");
        println!("Finished sw6");
    }
}
