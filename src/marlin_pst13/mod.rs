use crate::kzg10;
use crate::{BTreeMap, BTreeSet, MVPolynomial, String, Term, Vec};
use crate::{BatchLCProof, Error, Evaluations, QuerySet};
use crate::{LabeledCommitment, LabeledPolynomial, LinearCombination};
use crate::{PCRandomness, PCUniversalParams, PolynomialCommitment};
use algebra_core::msm::{FixedBaseMSM, VariableBaseMSM};
use algebra_core::{
    AffineCurve, One, PairingEngine, PrimeField, ProjectiveCurve, UniformRand, Zero,
};
use combinations::Combinations;
use core::{convert::TryInto, marker::PhantomData, ops::Index};
use rand_core::RngCore;

mod data_structures;
pub use data_structures::*;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Multivariate polynomial commitment based on the construction in [[PST13]][pst]
/// with batching and (optional) hiding property inspired by the univariate scheme
/// in [[CHMMVW20, "Marlin"]][marlin]
///
/// [pst]: https://eprint.iacr.org/2011/587.pdf
/// [marlin]: https://eprint.iacr.org/2019/104
pub struct MarlinPST13<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Domain: Index<usize, Output = E::Fr>,
{
    _engine: PhantomData<E>,
    _poly: PhantomData<P>,
}

impl<E, P> MarlinPST13<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Domain: Index<usize, Output = E::Fr>,
{
    /// MSM for `commitments` and `coeffs`
    fn combine_commitments<'a>(
        coeffs_and_comms: impl IntoIterator<Item = (E::Fr, &'a kzg10::Commitment<E>)>,
    ) -> E::G1Projective {
        let mut combined_comm = E::G1Projective::zero();
        for (coeff, comm) in coeffs_and_comms {
            if coeff.is_one() {
                combined_comm.add_assign_mixed(&comm.0);
            } else {
                combined_comm += &comm.0.mul(coeff);
            }
        }
        combined_comm
    }

    fn normalize_commitments<'a>(commitments: Vec<E::G1Projective>) -> Vec<kzg10::Commitment<E>> {
        let norm_comms = E::G1Projective::batch_normalization_into_affine(&commitments);
        norm_comms.into_iter().map(kzg10::Commitment).collect()
    }

    /// Accumulate `commitments` and `values` according to `opening_challenge`.
    fn accumulate_commitments_and_values<'a>(
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<kzg10::Commitment<E>>>,
        values: impl IntoIterator<Item = E::Fr>,
        opening_challenge: E::Fr,
    ) -> Result<(E::G1Projective, E::Fr), Error> {
        let acc_time = start_timer!(|| "Accumulating commitments and values");
        let mut combined_comm = E::G1Projective::zero();
        let mut combined_value = E::Fr::zero();
        let mut challenge_i = E::Fr::one();
        for (labeled_commitment, value) in commitments.into_iter().zip(values) {
            let commitment = labeled_commitment.commitment();
            challenge_i *= &opening_challenge;
            combined_comm += &commitment.0.mul(challenge_i);
            combined_value += &(value * &challenge_i);
        }
        end_timer!(acc_time);
        Ok((combined_comm, combined_value))
    }

    /// Check that the powers support the hiding bound
    fn check_hiding_bound(hiding_poly_degree: usize, num_powers: usize) -> Result<(), Error> {
        if hiding_poly_degree == 0 {
            Err(Error::HidingBoundIsZero)
        } else if hiding_poly_degree >= num_powers {
            // The above check uses `>=` because committing to a hiding poly with
            // degree `hiding_poly_degree` requires `hiding_poly_degree + 1`
            // powers.
            Err(Error::HidingBoundToolarge {
                hiding_poly_degree,
                num_powers,
            })
        } else {
            Ok(())
        }
    }

    /// Check that a given polynomial is supported by parameters
    fn check_degrees_and_bounds<'a>(
        supported_degree: usize,
        p: &'a LabeledPolynomial<'a, E::Fr, P>,
    ) -> Result<(), Error>
    where
        P: 'a,
    {
        if p.degree() > supported_degree {
            return Err(Error::PolynomialDegreeTooLarge {
                poly_degree: p.degree(),
                supported_degree,
                label: p.label().to_string(),
            });
        } else {
            Ok(())
        }
    }

    /// Convert polynomial coefficients to `BigInt`
    fn convert_to_bigints(p: &P) -> Vec<<E::Fr as PrimeField>::BigInt> {
        let to_bigint_time = start_timer!(|| "Converting polynomial coeffs to bigints");
        let plain_coeffs = ff_fft::cfg_into_iter!(p.terms())
            .map(|(_, coeff)| coeff.into_repr())
            .collect();
        end_timer!(to_bigint_time);
        plain_coeffs
    }
}

impl<E, P> PolynomialCommitment<E::Fr, P> for MarlinPST13<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr>,
    P::Domain: Index<usize, Output = E::Fr>,
{
    type UniversalParams = UniversalParams<E, P>;
    type CommitterKey = CommitterKey<E, P>;
    type VerifierKey = VerifierKey<E>;
    type Commitment = kzg10::Commitment<E>;
    type Randomness = Randomness<E, P>;
    type Proof = Proof<E>;
    type BatchProof = Vec<Self::Proof>;
    type Error = Error;

    /// Constructs public parameters when given as input the maximum degree `max_degree`
    /// and number of variables `num_vars` for the polynomial commitment scheme.
    fn setup<R: RngCore>(
        max_degree: usize,
        num_vars: Option<usize>,
        rng: &mut R,
    ) -> Result<UniversalParams<E, P>, Error> {
        let num_vars = num_vars.ok_or(Error::InvalidNumberOfVariables)?;
        if num_vars < 1 {
            return Err(Error::InvalidNumberOfVariables);
        }
        if max_degree < 1 {
            return Err(Error::DegreeIsZero);
        }
        let setup_time = start_timer!(|| format!(
            "MarlinPST13::Setup with {} variables and max degree {}",
            num_vars, max_degree
        ));
        // Trapdoor evaluation points
        let mut betas = Vec::with_capacity(num_vars);
        for _ in 0..num_vars {
            betas.push(E::Fr::rand(rng));
        }
        // Generators
        let g = E::G1Projective::rand(rng);
        let gamma_g = E::G1Projective::rand(rng);
        let h = E::G2Projective::rand(rng);

        // A list of all variable numbers of multiplicity `max_degree`
        let variable_set: Vec<_> = (0..num_vars)
            .flat_map(|var| vec![var; max_degree])
            .collect();
        // Generate all possible monomials with `1 <= degree <= max_degree`
        let (powers_of_beta, mut powers_of_beta_terms): (Vec<_>, Vec<_>) =
            ff_fft::cfg_into_iter!(1..=max_degree)
                .flat_map(|degree| {
                    // Sample all combinations of `degree` variables from `variable_set`
                    let terms: Vec<Vec<usize>> = if variable_set.len() == degree {
                        vec![variable_set.clone()]
                    } else {
                        Combinations::new(variable_set.clone(), degree).collect()
                    };
                    // For each multiset in `terms` evaluate the corresponding monomial at the
                    // trapdoor and generate a `P::Term` object to index it
                    terms
                        .into_iter()
                        .map(|term| {
                            let value: E::Fr = term.iter().map(|e| betas[*e]).product();
                            let term = (0..num_vars)
                                .map(|var| (var, term.iter().filter(|e| **e == var).count()))
                                .collect();
                            (value, P::Term::new(term))
                        })
                        .collect::<Vec<_>>()
                })
                .unzip();

        let scalar_bits = E::Fr::size_in_bits();
        let g_time = start_timer!(|| "Generating powers of G");
        let window_size = FixedBaseMSM::get_mul_window_size(max_degree + 1);
        let g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, g);
        let mut powers_of_g = FixedBaseMSM::multi_scalar_mul::<E::G1Projective>(
            scalar_bits,
            window_size,
            &g_table,
            &powers_of_beta,
        );
        powers_of_g.push(g);
        powers_of_beta_terms.push(P::Term::new(vec![]));
        end_timer!(g_time);

        let gamma_g_time = start_timer!(|| "Generating powers of gamma * G");
        let window_size = FixedBaseMSM::get_mul_window_size(max_degree + 2);
        let gamma_g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, gamma_g);
        // Each element `i` of `powers_of_gamma_g` is a vector of length `max_degree+1`
        // containing `betas[i]^j \gamma G` for `j` from 1 to `max_degree+1` to support
        // up to `max_degree` queries
        let mut powers_of_gamma_g = vec![Vec::new(); num_vars];
        ff_fft::cfg_iter_mut!(powers_of_gamma_g)
            .enumerate()
            .for_each(|(i, v)| {
                let mut powers_of_beta = Vec::with_capacity(max_degree);
                let mut cur = E::Fr::one();
                for _ in 0..=max_degree {
                    cur *= &betas[i];
                    powers_of_beta.push(cur);
                }
                *v = FixedBaseMSM::multi_scalar_mul::<E::G1Projective>(
                    scalar_bits,
                    window_size,
                    &gamma_g_table,
                    &powers_of_beta,
                );
            });
        end_timer!(gamma_g_time);

        let powers_of_g = E::G1Projective::batch_normalization_into_affine(&powers_of_g);
        let gamma_g = gamma_g.into_affine();
        let powers_of_gamma_g = powers_of_gamma_g
            .into_iter()
            .map(|v| E::G1Projective::batch_normalization_into_affine(&v))
            .collect();
        let beta_h: Vec<_> = betas.iter().map(|b| h.mul(*b).into_affine()).collect();
        let h = h.into_affine();
        let prepared_h = h.into();
        let prepared_beta_h = beta_h.iter().map(|bh| (*bh).into()).collect();

        // Convert `powers_of_g` to a BTreeMap indexed by `powers_of_beta_terms`
        let powers_of_g = powers_of_beta_terms
            .into_iter()
            .zip(powers_of_g.into_iter())
            .collect();

        let pp = UniversalParams {
            num_vars,
            max_degree,
            powers_of_g,
            gamma_g,
            powers_of_gamma_g,
            h,
            beta_h,
            prepared_h,
            prepared_beta_h,
        };
        end_timer!(setup_time);
        Ok(pp)
    }

    /// Specializes the public parameters for polynomials up to the given `supported_degree`
    ///
    /// TODO: Add the ability to trim the number of variables
    /// TODO: Update for support_hiding_bound
    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        _supported_hiding_bound: usize,
        _enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        let max_degree = pp.max_degree();
        if supported_degree > max_degree {
            return Err(Error::TrimmingDegreeTooLarge);
        }

        let ck_time = start_timer!(|| format!(
            "Constructing CommitterKey of size {} for unshifted polys",
            supported_degree
        ));
        // We want to support making up to supported_degree queries to committed
        // polynomials.
        let powers_of_g = pp
            .powers_of_g
            .iter()
            .filter(|(k, _)| k.degree() <= supported_degree)
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        let powers_of_gamma_g = pp
            .powers_of_gamma_g
            .iter()
            .map(|e| e[..=supported_degree].to_vec())
            .collect();
        end_timer!(ck_time);

        let ck = CommitterKey {
            powers_of_g,
            gamma_g: pp.gamma_g,
            powers_of_gamma_g,
            num_vars: pp.num_vars,
            supported_degree,
            max_degree,
        };

        let vk = VerifierKey {
            g: pp.powers_of_g[&P::Term::new(vec![])],
            gamma_g: pp.gamma_g,
            h: pp.h,
            beta_h: pp.beta_h.clone(),
            prepared_h: pp.prepared_h.clone(),
            prepared_beta_h: pp.prepared_beta_h.clone(),
            num_vars: pp.num_vars,
            supported_degree,
            max_degree,
        };
        Ok((ck, vk))
    }

    /// Outputs a commitments to `polynomials`.
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, E::Fr, P>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    >
    where
        P: 'a,
    {
        let rng = &mut crate::optional_rng::OptionalRng(rng);
        let commit_time = start_timer!(|| "Committing to polynomials");
        let mut commitments = Vec::new();
        let mut randomness = Vec::new();
        for p in polynomials {
            let label = p.label();
            let hiding_bound = p.hiding_bound();
            let polynomial: &P = p.polynomial();
            Self::check_degrees_and_bounds(ck.supported_degree, &p)?;

            let commit_time = start_timer!(|| {
                format!(
                    "Polynomial {} with {} variables of degree {} and hiding bound {:?}",
                    label,
                    polynomial.num_vars,
                    polynomial.degree(),
                    hiding_bound,
                )
            });
            // Get the powers of `G` corresponding to the terms of `polynomial`
            let powers_of_g = polynomial
                .terms()
                .iter()
                .map(|(term, _)| *ck.powers_of_g.get(term).unwrap())
                .collect::<Vec<_>>();
            // Convert coefficients of `polynomial` to BigInts
            let plain_ints = Self::convert_to_bigints(&polynomial);

            let msm_time = start_timer!(|| "MSM to compute commitment to plaintext poly");
            let mut commitment = VariableBaseMSM::multi_scalar_mul(&powers_of_g, &plain_ints);
            end_timer!(msm_time);

            // Sample random polynomial
            let mut rand = Randomness::<E, P>::empty();
            if let Some(hiding_degree) = hiding_bound {
                let sample_random_poly_time = start_timer!(|| format!(
                    "Sampling a random polynomial of degree {}",
                    hiding_degree
                ));
                rand = Randomness::rand(hiding_degree, false, Some(ck.num_vars), rng);
                Self::check_hiding_bound(hiding_degree, ck.supported_degree + 1)?;
                end_timer!(sample_random_poly_time);
            }

            // Get the powers of `\gamma G` corresponding to the terms of `rand`
            let powers_of_gamma_g = rand
                .blinding_polynomial
                .terms()
                .iter()
                .map(|(v, _)| {
                    // Implicit Assumption: Each monomial in `rand` is univariate
                    let vars = v.vars();
                    match v.is_constant() {
                        true => ck.gamma_g,
                        false => ck.powers_of_gamma_g[vars[0]][v.degree() - 1],
                    }
                })
                .collect::<Vec<_>>();
            // Convert coefficients of `rand` to BigInt
            let random_ints = Self::convert_to_bigints(&rand.blinding_polynomial);

            let msm_time = start_timer!(|| "MSM to compute commitment to random poly");
            let random_commitment =
                VariableBaseMSM::multi_scalar_mul(&powers_of_gamma_g, &random_ints).into_affine();
            end_timer!(msm_time);

            // Mask commitment with random poly
            commitment.add_assign_mixed(&random_commitment);

            let comm = kzg10::Commitment(commitment.into());
            commitments.push(LabeledCommitment::new(label.to_string(), comm, None));
            randomness.push(rand);
            end_timer!(commit_time);
        }
        end_timer!(commit_time);
        Ok((commitments, randomness))
    }

    /// Output a proof of evaluation of `labeled_polynomials` at `point`
    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, E::Fr, P>>,
        _commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &P::Domain,
        opening_challenge: E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        Self::Randomness: 'a,
        Self::Commitment: 'a,
        P: 'a,
    {
        // Compute random linear combinations of committed polynomials and randomness
        let mut p = P::zero();
        let mut r = Randomness::empty();
        let mut challenge_j = E::Fr::one();
        for (polynomial, rand) in labeled_polynomials.into_iter().zip(rands) {
            Self::check_degrees_and_bounds(ck.supported_degree, &polynomial)?;
            challenge_j *= &opening_challenge;
            p += (challenge_j, polynomial.polynomial());
            r += (challenge_j, rand);
        }

        let open_time = start_timer!(|| format!("Opening polynomial of degree {}", p.degree()));
        let witness_time = start_timer!(|| "Computing witness polynomials");
        let witnesses = p.divide_at_point(point);
        let hiding_witnesses = if r.is_hiding() {
            Some(r.blinding_polynomial.divide_at_point(point))
        } else {
            None
        };
        end_timer!(witness_time);

        let witness_comm_time = start_timer!(|| "Computing commitment to witness polynomials");
        let mut w = ff_fft::cfg_into_iter!(witnesses)
            .map(|w| {
                // Get the powers of `G` corresponding to the witness poly
                let powers_of_g = w
                    .terms()
                    .iter()
                    .map(|(term, _)| *ck.powers_of_g.get(term).unwrap())
                    .collect::<Vec<_>>();
                // Convert coefficients to BigInt
                let witness_ints = Self::convert_to_bigints(&w);
                // Compute MSM
                VariableBaseMSM::multi_scalar_mul(&powers_of_g, &witness_ints)
            })
            .collect::<Vec<_>>();
        end_timer!(witness_comm_time);

        // If the evaluation should be hiding, compute the MSM for `hiding_witnesses` and add
        // to the `w`. Additionally, compute the evaluation of `r` at `point`.
        let random_v = if let Some(hiding_witnesses) = hiding_witnesses {
            let witness_comm_time =
                start_timer!(|| "Computing commitment to hiding witness polynomials");
            ff_fft::cfg_iter_mut!(w)
                .enumerate()
                .for_each(|(i, witness)| {
                    let hiding_witness = &hiding_witnesses[i];
                    // Get the powers of `\gamma G` corresponding to the terms of `hiding_witness`
                    let powers_of_gamma_g = hiding_witness
                        .terms()
                        .iter()
                        .map(|(v, _)| {
                            // Implicit Assumption: Each monomial in `hiding_witness` is univariate
                            let vars = v.vars();
                            match v.is_constant() {
                                true => ck.gamma_g,
                                false => ck.powers_of_gamma_g[vars[0]][v.degree() - 1],
                            }
                        })
                        .collect::<Vec<_>>();
                    // Convert coefficients to BigInt
                    let hiding_witness_ints = Self::convert_to_bigints(hiding_witness);
                    // Compute MSM and add result to witness
                    *witness += &VariableBaseMSM::multi_scalar_mul(
                        &powers_of_gamma_g,
                        &hiding_witness_ints,
                    );
                });
            end_timer!(witness_comm_time);
            Some(r.blinding_polynomial.evaluate(point))
        } else {
            None
        };
        end_timer!(open_time);
        Ok(Proof {
            w: w.into_iter().map(|w| w.into_affine()).collect(),
            random_v,
        })
    }

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    fn check<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &P::Domain,
        values: impl IntoIterator<Item = E::Fr>,
        proof: &Self::Proof,
        opening_challenge: E::Fr,
        _rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        // TODO: More comments to benchmark individual speeds here
        let check_time = start_timer!(|| "Checking evaluations");
        // Accumulate commitments and values
        let (combined_comm, combined_value) =
            Self::accumulate_commitments_and_values(commitments, values, opening_challenge)?;
        // Compute both sides of the pairing equation
        let mut inner = combined_comm.into().into_projective() - &vk.g.mul(combined_value);
        if let Some(random_v) = proof.random_v {
            inner -= &vk.gamma_g.mul(random_v);
        }
        let lhs = E::pairing(inner, vk.h);

        // Create a list of elements corresponding to each pairing in the product on the rhs
        let rhs_product: Vec<(E::G1Prepared, E::G2Prepared)> = proof
            .w
            .iter()
            .enumerate()
            .map(|(j, w_j)| {
                let beta_minus_z: E::G2Affine =
                    (vk.beta_h[j].into_projective() - &vk.h.mul(point[j])).into();
                ((*w_j).into(), beta_minus_z.into())
            })
            .collect();
        let rhs = E::product_of_pairings(&rhs_product);
        end_timer!(check_time);

        Ok(lhs == rhs)
    }

    fn batch_check<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Domain>,
        values: &Evaluations<E::Fr, P::Domain>,
        proof: &Self::BatchProof,
        opening_challenge: E::Fr,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let commitments: BTreeMap<_, _> = commitments.into_iter().map(|c| (c.label(), c)).collect();
        let mut query_to_labels_map = BTreeMap::new();

        for (label, point) in query_set.iter() {
            let labels = query_to_labels_map.entry(point).or_insert(BTreeSet::new());
            labels.insert(label);
        }
        assert_eq!(proof.len(), query_to_labels_map.len());

        let mut combined_comms = Vec::new();
        let mut combined_queries = Vec::new();
        let mut combined_evals = Vec::new();
        for (query, labels) in query_to_labels_map.into_iter() {
            let lc_time =
                start_timer!(|| format!("Randomly combining {} commitments", labels.len()));
            let mut comms_to_combine: Vec<&'_ LabeledCommitment<_>> = Vec::new();
            let mut values_to_combine = Vec::new();
            for label in labels.into_iter() {
                let commitment = commitments.get(label).ok_or(Error::MissingPolynomial {
                    label: label.to_string(),
                })?;

                let v_i = values.get(&(label.clone(), query.clone())).ok_or(
                    Error::MissingEvaluation {
                        label: label.to_string(),
                    },
                )?;

                comms_to_combine.push(commitment);
                values_to_combine.push(*v_i);
            }
            let (c, v) = Self::accumulate_commitments_and_values(
                comms_to_combine,
                values_to_combine,
                opening_challenge,
            )?;
            end_timer!(lc_time);
            combined_comms.push(c);
            combined_queries.push(query);
            combined_evals.push(v);
        }
        let norm_time = start_timer!(|| "Normalizaing combined commitments");
        E::G1Projective::batch_normalization(&mut combined_comms);
        let combined_comms = combined_comms
            .into_iter()
            .map(|c| kzg10::Commitment::<E>(c.into()))
            .collect::<Vec<_>>();
        end_timer!(norm_time);

        let check_time =
            start_timer!(|| format!("Checking {} evaluation proofs", commitments.len()));
        let g = vk.g.into_projective();
        let gamma_g = vk.gamma_g.into_projective();

        let mut total_c = <E::G1Projective>::zero();
        let mut total_w = vec![<E::G1Projective>::zero(); vk.num_vars];

        let combination_time = start_timer!(|| "Combining commitments and proofs");
        let mut randomizer = E::Fr::one();
        // Instead of multiplying g and gamma_g in each turn, we simply accumulate
        // their coefficients and perform a final multiplication at the end.
        let mut g_multiplier = E::Fr::zero();
        let mut gamma_g_multiplier = E::Fr::zero();
        for (((c, z), v), proof) in combined_comms
            .iter()
            .zip(combined_queries)
            .zip(combined_evals)
            .zip(proof)
        {
            let w = &proof.w;
            let mut temp: E::G1Projective =
                w.iter().enumerate().map(|(j, w_j)| w_j.mul(z[j])).sum();
            temp.add_assign_mixed(&c.0);
            let c = temp;
            g_multiplier += &(randomizer * &v);
            if let Some(random_v) = proof.random_v {
                gamma_g_multiplier += &(randomizer * &random_v);
            }
            total_c += &c.mul(randomizer);
            total_w
                .iter_mut()
                .enumerate()
                .for_each(|(i, w_i)| *w_i += &w[i].mul(randomizer));
            // We don't need to sample randomizers from the full field,
            // only from 128-bit strings.
            randomizer = u128::rand(rng).into();
        }
        total_c -= &g.mul(g_multiplier);
        total_c -= &gamma_g.mul(gamma_g_multiplier);
        end_timer!(combination_time);

        // Converting results for pairing
        let mut pairings = Vec::new();
        total_w.into_iter().enumerate().for_each(|(j, w_j)| {
            pairings.push(((-w_j).into_affine().into(), vk.prepared_beta_h[j].clone()))
        });
        pairings.push((total_c.into_affine().into(), vk.prepared_h.clone()));

        let pairing_time = start_timer!(|| "Performing product of pairings");
        let result = E::product_of_pairings(&pairings).is_one();
        end_timer!(pairing_time);
        end_timer!(check_time);
        Ok(result)
    }

    fn open_combinations<'a>(
        ck: &Self::CommitterKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<E::Fr>>,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<'a, E::Fr, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Domain>,
        opening_challenge: E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<BatchLCProof<E::Fr, P, Self>, Self::Error>
    where
        Self::Randomness: 'a,
        Self::Commitment: 'a,
        P: 'a,
    {
        let label_map = polynomials
            .into_iter()
            .zip(rands)
            .zip(commitments)
            .map(|((p, r), c)| (p.label(), (p, r, c)))
            .collect::<BTreeMap<_, _>>();

        let mut lc_polynomials = Vec::new();
        let mut lc_randomness = Vec::new();
        let mut lc_commitments = Vec::new();
        let mut lc_info = Vec::new();

        for lc in lc_s {
            let lc_label = lc.label().clone();
            let mut poly = P::zero();
            let mut hiding_bound = None;
            let mut randomness = Self::Randomness::empty();
            let mut coeffs_and_comms = Vec::new();

            for (coeff, label) in lc.iter().filter(|(_, l)| !l.is_one()) {
                let label: &String = label.try_into().expect("cannot be one!");
                let &(cur_poly, cur_rand, cur_comm) =
                    label_map.get(label).ok_or(Error::MissingPolynomial {
                        label: label.to_string(),
                    })?;

                // Some(_) > None, always.
                hiding_bound = core::cmp::max(hiding_bound, cur_poly.hiding_bound());
                poly += (*coeff, cur_poly.polynomial());
                randomness += (*coeff, cur_rand);
                coeffs_and_comms.push((*coeff, cur_comm.commitment()));
            }
            let lc_poly = LabeledPolynomial::new_owned(lc_label.clone(), poly, None, hiding_bound);
            lc_polynomials.push(lc_poly);
            lc_randomness.push(randomness);
            lc_commitments.push(Self::combine_commitments(coeffs_and_comms));
            lc_info.push((lc_label, None));
        }

        let comms = Self::normalize_commitments(lc_commitments);
        let lc_commitments = lc_info
            .into_iter()
            .zip(comms)
            .map(|((label, d), c)| LabeledCommitment::new(label, c, d))
            .collect::<Vec<_>>();

        let proof = Self::batch_open(
            ck,
            lc_polynomials.iter(),
            lc_commitments.iter(),
            &query_set,
            opening_challenge,
            lc_randomness.iter(),
            rng,
        )?;

        Ok(BatchLCProof { proof, evals: None })
    }

    /// Checks that `values` are the true evaluations at `query_set` of the polynomials
    /// committed in `labeled_commitments`.
    fn check_combinations<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<E::Fr>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Domain>,
        evaluations: &Evaluations<E::Fr, P::Domain>,
        proof: &BatchLCProof<E::Fr, P, Self>,
        opening_challenge: E::Fr,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let BatchLCProof { proof, .. } = proof;
        let label_comm_map = commitments
            .into_iter()
            .map(|c| (c.label(), c))
            .collect::<BTreeMap<_, _>>();

        let mut lc_commitments = Vec::new();
        let mut lc_info = Vec::new();
        let mut evaluations = evaluations.clone();

        let lc_processing_time = start_timer!(|| "Combining commitments");
        for lc in lc_s {
            let lc_label = lc.label().clone();
            let mut coeffs_and_comms = Vec::new();
            for (coeff, label) in lc.iter() {
                if label.is_one() {
                    for (&(ref label, _), ref mut eval) in evaluations.iter_mut() {
                        if label == &lc_label {
                            **eval -= coeff;
                        }
                    }
                } else {
                    let label: &String = label.try_into().unwrap();
                    let &cur_comm = label_comm_map.get(label).ok_or(Error::MissingPolynomial {
                        label: label.to_string(),
                    })?;
                    coeffs_and_comms.push((*coeff, cur_comm.commitment()));
                }
            }
            let lc_time = start_timer!(|| format!("Combining commitments for {}", lc_label));
            lc_commitments.push(Self::combine_commitments(coeffs_and_comms));
            end_timer!(lc_time);
            lc_info.push((lc_label, None));
        }
        end_timer!(lc_processing_time);
        let combined_comms_norm_time = start_timer!(|| "Normalizing commitments");
        let comms = Self::normalize_commitments(lc_commitments);
        let lc_commitments = lc_info
            .into_iter()
            .zip(comms)
            .map(|((label, d), c)| LabeledCommitment::new(label, c, d))
            .collect::<Vec<_>>();
        end_timer!(combined_comms_norm_time);

        Self::batch_check(
            vk,
            &lc_commitments,
            &query_set,
            &evaluations,
            proof,
            opening_challenge,
            rng,
        )
    }
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]
    use super::MarlinPST13;
    use algebra::Bls12_377;
    use algebra::Bls12_381;
    use algebra::PairingEngine;

    use ff_fft::multivariate::{SparsePolynomial as SparsePoly, SparseTerm};

    type MVPoly_381 = SparsePoly<<Bls12_381 as PairingEngine>::Fr, SparseTerm>;
    type MVPoly_377 = SparsePoly<<Bls12_377 as PairingEngine>::Fr, SparseTerm>;

    type PC<E, P> = MarlinPST13<E, P>;
    type PC_Bls12_381 = PC<Bls12_381, MVPoly_381>;
    type PC_Bls12_377 = PC<Bls12_377, MVPoly_377>;

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        let num_vars = Some(10);
        single_poly_test::<_, _, PC_Bls12_377>(num_vars).expect("test failed for bls12-377");
        single_poly_test::<_, _, PC_Bls12_381>(num_vars).expect("test failed for bls12-381");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        let num_vars = Some(10);
        full_end_to_end_test::<_, _, PC_Bls12_377>(num_vars).expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_test::<_, _, PC_Bls12_381>(num_vars).expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        let num_vars = Some(10);
        single_equation_test::<_, _, PC_Bls12_377>(num_vars).expect("test failed for bls12-377");
        println!("Finished bls12-377");
        single_equation_test::<_, _, PC_Bls12_381>(num_vars).expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        let num_vars = Some(10);
        two_equation_test::<_, _, PC_Bls12_377>(num_vars).expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_test::<_, _, PC_Bls12_381>(num_vars).expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        let num_vars = Some(10);
        full_end_to_end_equation_test::<_, _, PC_Bls12_377>(num_vars)
            .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_equation_test::<_, _, PC_Bls12_381>(num_vars)
            .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }
}
