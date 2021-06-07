use crate::{
    kzg10,
    marlin::{marlin_pc, Marlin},
};
use crate::{BatchLCProof, Error, Evaluations, QuerySet};
use crate::{LabeledCommitment, LabeledPolynomial, LinearCombination};
use crate::{PCRandomness, PCUniversalParams, PolynomialCommitment};
use crate::{ToString, Vec};
use ark_ec::{
    msm::{FixedBaseMSM, VariableBaseMSM},
    AffineCurve, PairingEngine, ProjectiveCurve,
};
use ark_ff::{One, PrimeField, UniformRand, Zero};
use ark_poly::{multivariate::Term, MVPolynomial};
use ark_std::rand::RngCore;
use ark_std::{marker::PhantomData, ops::Index, vec};

mod data_structures;
pub use data_structures::*;

mod combinations;
use combinations::*;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Multivariate polynomial commitment based on the construction in [[PST13]][pst]
/// with batching and (optional) hiding property inspired by the univariate scheme
/// in [[CHMMVW20, "Marlin"]][marlin]
///
/// [pst]: https://eprint.iacr.org/2011/587
/// [marlin]: https://eprint.iacr.org/2019/104
pub struct MarlinPST13<E: PairingEngine, P: MVPolynomial<E::Fr>> {
    _engine: PhantomData<E>,
    _poly: PhantomData<P>,
}

impl<E: PairingEngine, P: MVPolynomial<E::Fr>> MarlinPST13<E, P> {
    /// Given some point `z`, compute the quotients `w_i(X)` s.t
    ///
    /// `p(X) - p(z) = (X_1-z_1)*w_1(X) + (X_2-z_2)*w_2(X) + ... + (X_l-z_l)*w_l(X)`
    ///
    /// These quotients can always be found with no remainder.
    fn divide_at_point(p: &P, point: &P::Point) -> Vec<P>
    where
        P::Point: Index<usize, Output = E::Fr>,
    {
        let num_vars = p.num_vars();
        if p.is_zero() {
            return vec![P::zero(); num_vars];
        }
        let mut quotients = Vec::with_capacity(num_vars);
        // `cur` represents the current dividend
        let mut cur = p.clone();
        // Divide `cur` by `X_i - z_i`
        for i in 0..num_vars {
            let mut quotient_terms = Vec::new();
            let mut remainder_terms = Vec::new();
            for (mut coeff, term) in cur.terms() {
                // Since the final remainder is guaranteed to be 0, all the constant terms
                // cancel out so we don't need to keep track of them
                if term.is_constant() {
                    continue;
                }
                // If the current term contains `X_i` then divide appropiately,
                // otherwise add it to the remainder
                let mut term_vec = (&*term).to_vec();
                match term_vec.binary_search_by(|(var, _)| var.cmp(&i)) {
                    Ok(idx) => {
                        // Repeatedly divide the term by `X_i - z_i` until the remainder
                        // doesn't contain any `X_i`s
                        while term_vec[idx].1 > 1 {
                            // First divide by `X_i` and add the term to the quotient
                            term_vec[idx] = (i, term_vec[idx].1 - 1);
                            quotient_terms.push((coeff, P::Term::new(term_vec.clone())));
                            // Then compute the remainder term in-place
                            coeff *= &point[i];
                        }
                        // Since `X_i` is power 1, we can remove it entirely
                        term_vec.remove(idx);
                        quotient_terms.push((coeff, P::Term::new(term_vec.clone())));
                        remainder_terms.push((point[i] * &coeff, P::Term::new(term_vec)));
                    }
                    Err(_) => remainder_terms.push((coeff, term.clone())),
                }
            }
            quotients.push(P::from_coefficients_vec(num_vars, quotient_terms));
            // Set the current dividend to be the remainder of this division
            cur = P::from_coefficients_vec(num_vars, remainder_terms);
        }
        quotients
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
        p: &'a LabeledPolynomial<E::Fr, P>,
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
        let plain_coeffs = ark_std::cfg_into_iter!(p.terms())
            .map(|(coeff, _)| coeff.into_repr())
            .collect();
        plain_coeffs
    }
}

impl<E, P> PolynomialCommitment<E::Fr, P> for MarlinPST13<E, P>
where
    E: PairingEngine,
    P: MVPolynomial<E::Fr> + Sync,
    P::Point: Index<usize, Output = E::Fr>,
{
    type UniversalParams = UniversalParams<E, P>;
    type CommitterKey = CommitterKey<E, P>;
    type VerifierKey = VerifierKey<E>;
    type PreparedVerifierKey = PreparedVerifierKey<E>;
    type Commitment = marlin_pc::Commitment<E>;
    type PreparedCommitment = marlin_pc::PreparedCommitment<E>;
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
        let (powers_of_beta, mut powers_of_beta_terms): (Vec<_>, Vec<_>) = (1..=max_degree)
            .flat_map(|degree| {
                // Sample all combinations of `degree` variables from `variable_set`
                let terms: Vec<Vec<usize>> = if variable_set.len() == degree {
                    vec![variable_set.clone()]
                } else {
                    Combinations::new(variable_set.clone(), degree).collect()
                };
                // For each multiset in `terms` evaluate the corresponding monomial at the
                // trapdoor and generate a `P::Term` object to index it
                ark_std::cfg_into_iter!(terms)
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
        ark_std::cfg_iter_mut!(powers_of_gamma_g)
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
        let beta_h: Vec<_> = betas
            .iter()
            .map(|b| h.mul(&(*b).into_repr()).into_affine())
            .collect();
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
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::Fr, P>>,
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
                    "Polynomial {} with degree {} and hiding bound {:?}",
                    label,
                    polynomial.degree(),
                    hiding_bound,
                )
            });
            // Get the powers of `G` corresponding to the terms of `polynomial`
            let powers_of_g = ark_std::cfg_iter!(polynomial.terms())
                .map(|(_, term)| *ck.powers_of_g.get(term).unwrap())
                .collect::<Vec<_>>();
            // Convert coefficients of `polynomial` to BigInts
            let to_bigint_time = start_timer!(|| "Converting polynomial coeffs to bigints");
            let plain_ints = Self::convert_to_bigints(&polynomial);
            end_timer!(to_bigint_time);

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
                .map(|(_, term)| {
                    // Implicit Assumption: Each monomial in `rand` is univariate
                    let vars = term.vars();
                    match term.is_constant() {
                        true => ck.gamma_g,
                        false => ck.powers_of_gamma_g[vars[0]][term.degree() - 1],
                    }
                })
                .collect::<Vec<_>>();
            // Convert coefficients of `rand` to BigInt
            let to_bigint_time = start_timer!(|| "Converting polynomial coeffs to bigints");
            let random_ints = Self::convert_to_bigints(&rand.blinding_polynomial);
            end_timer!(to_bigint_time);

            let msm_time = start_timer!(|| "MSM to compute commitment to random poly");
            let random_commitment =
                VariableBaseMSM::multi_scalar_mul(&powers_of_gamma_g, &random_ints).into_affine();
            end_timer!(msm_time);

            // Mask commitment with random poly
            commitment.add_assign_mixed(&random_commitment);

            let comm = Self::Commitment {
                comm: kzg10::Commitment(commitment.into()),
                shifted_comm: None,
            };

            commitments.push(LabeledCommitment::new(label.to_string(), comm, None));
            randomness.push(rand);
            end_timer!(commit_time);
        }
        end_timer!(commit_time);
        Ok((commitments, randomness))
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the same.
    fn open_individual_opening_challenges<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::Fr, P>>,
        _commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &P::Point,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        P: 'a,
        Self::Randomness: 'a,
        Self::Commitment: 'a,
    {
        // Compute random linear combinations of committed polynomials and randomness
        let mut p = P::zero();
        let mut r = Randomness::empty();
        let mut opening_challenge_counter = 0;
        for (polynomial, rand) in labeled_polynomials.into_iter().zip(rands) {
            Self::check_degrees_and_bounds(ck.supported_degree, &polynomial)?;

            // compute challenge^j and challenge^{j+1}.
            let challenge_j = opening_challenges(opening_challenge_counter);
            opening_challenge_counter += 1;

            p += (challenge_j, polynomial.polynomial());
            r += (challenge_j, rand);
        }

        let open_time = start_timer!(|| format!("Opening polynomial of degree {}", p.degree()));
        let witness_time = start_timer!(|| "Computing witness polynomials");
        let witnesses = Self::divide_at_point(&p, point);
        let hiding_witnesses = if r.is_hiding() {
            Some(Self::divide_at_point(&r.blinding_polynomial, point))
        } else {
            None
        };
        end_timer!(witness_time);

        let witness_comm_time = start_timer!(|| "Computing commitment to witness polynomials");
        let mut w = witnesses
            .iter()
            .map(|w| {
                // Get the powers of `G` corresponding to the witness poly
                let powers_of_g = ark_std::cfg_iter!(w.terms())
                    .map(|(_, term)| *ck.powers_of_g.get(term).unwrap())
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
            ark_std::cfg_iter_mut!(w)
                .enumerate()
                .for_each(|(i, witness)| {
                    let hiding_witness = &hiding_witnesses[i];
                    // Get the powers of `\gamma G` corresponding to the terms of `hiding_witness`
                    let powers_of_gamma_g = hiding_witness
                        .terms()
                        .iter()
                        .map(|(_, term)| {
                            // Implicit Assumption: Each monomial in `hiding_witness` is univariate
                            let vars = term.vars();
                            match term.is_constant() {
                                true => ck.gamma_g,
                                false => ck.powers_of_gamma_g[vars[0]][term.degree() - 1],
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
    fn check_individual_opening_challenges<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a P::Point,
        values: impl IntoIterator<Item = E::Fr>,
        proof: &Self::Proof,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let check_time = start_timer!(|| "Checking evaluations");
        // Accumulate commitments and values
        let (combined_comm, combined_value) =
            Marlin::accumulate_commitments_and_values_individual_opening_challenges(
                commitments,
                values,
                opening_challenges,
                None,
            )?;
        // Compute both sides of the pairing equation
        let mut inner = combined_comm.into().into_projective() - &vk.g.mul(combined_value);
        if let Some(random_v) = proof.random_v {
            inner -= &vk.gamma_g.mul(random_v);
        }
        let lhs = E::pairing(inner, vk.h);

        // Create a list of elements corresponding to each pairing in the product on the rhs
        let rhs_product: Vec<(E::G1Prepared, E::G2Prepared)> = ark_std::cfg_iter!(proof.w)
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

    fn batch_check_individual_opening_challenges<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        values: &Evaluations<P::Point, E::Fr>,
        proof: &Self::BatchProof,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let (combined_comms, combined_queries, combined_evals) = Marlin::combine_and_normalize(
            commitments,
            query_set,
            values,
            opening_challenges,
            None,
        )?;
        let check_time =
            start_timer!(|| format!("Checking {} evaluation proofs", combined_comms.len()));
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
            let mut temp: E::G1Projective = ark_std::cfg_iter!(w)
                .enumerate()
                .map(|(j, w_j)| w_j.mul(z[j]))
                .sum();
            temp.add_assign_mixed(&c.0);
            let c = temp;
            g_multiplier += &(randomizer * &v);
            if let Some(random_v) = proof.random_v {
                gamma_g_multiplier += &(randomizer * &random_v);
            }
            total_c += &c.mul(&randomizer.into_repr());
            ark_std::cfg_iter_mut!(total_w)
                .enumerate()
                .for_each(|(i, w_i)| *w_i += &w[i].mul(randomizer));
            // We don't need to sample randomizers from the full field,
            // only from 128-bit strings.
            randomizer = u128::rand(rng).into();
        }
        total_c -= &g.mul(&g_multiplier.into_repr());
        total_c -= &gamma_g.mul(&gamma_g_multiplier.into_repr());
        end_timer!(combination_time);

        let to_affine_time = start_timer!(|| "Converting results to affine for pairing");
        let mut pairings = Vec::new();
        total_w.into_iter().enumerate().for_each(|(j, w_j)| {
            pairings.push(((-w_j).into_affine().into(), vk.prepared_beta_h[j].clone()))
        });
        pairings.push((total_c.into_affine().into(), vk.prepared_h.clone()));
        end_timer!(to_affine_time);

        let pairing_time = start_timer!(|| "Performing product of pairings");
        let result = E::product_of_pairings(&pairings).is_one();
        end_timer!(pairing_time);
        end_timer!(check_time);
        Ok(result)
    }

    fn open_combinations_individual_opening_challenges<'a>(
        ck: &Self::CommitterKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<E::Fr>>,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<E::Fr, P>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<BatchLCProof<E::Fr, P, Self>, Self::Error>
    where
        P: 'a,
        Self::Randomness: 'a,
        Self::Commitment: 'a,
    {
        Marlin::open_combinations_individual_opening_challenges(
            ck,
            lc_s,
            polynomials,
            commitments,
            query_set,
            opening_challenges,
            rands,
            rng,
        )
    }

    /// Checks that `values` are the true evaluations at `query_set` of the polynomials
    /// committed in `labeled_commitments`.
    fn check_combinations_individual_opening_challenges<'a, R: RngCore>(
        vk: &Self::VerifierKey,
        lc_s: impl IntoIterator<Item = &'a LinearCombination<E::Fr>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        query_set: &QuerySet<P::Point>,
        evaluations: &Evaluations<P::Point, E::Fr>,
        proof: &BatchLCProof<E::Fr, P, Self>,
        opening_challenges: &dyn Fn(u64) -> E::Fr,
        rng: &mut R,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        Marlin::check_combinations_individual_opening_challenges(
            vk,
            lc_s,
            commitments,
            query_set,
            evaluations,
            proof,
            opening_challenges,
            rng,
        )
    }
}

#[cfg(test)]
mod tests {
    #![allow(non_camel_case_types)]
    use super::MarlinPST13;
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_ec::PairingEngine;
    use ark_ff::UniformRand;
    use ark_poly::{
        multivariate::{SparsePolynomial as SparsePoly, SparseTerm},
        MVPolynomial,
    };
    use ark_std::rand::rngs::StdRng;

    type MVPoly_381 = SparsePoly<<Bls12_381 as PairingEngine>::Fr, SparseTerm>;
    type MVPoly_377 = SparsePoly<<Bls12_377 as PairingEngine>::Fr, SparseTerm>;

    type PC<E, P> = MarlinPST13<E, P>;
    type PC_Bls12_381 = PC<Bls12_381, MVPoly_381>;
    type PC_Bls12_377 = PC<Bls12_377, MVPoly_377>;

    fn rand_poly<E: PairingEngine>(
        degree: usize,
        num_vars: Option<usize>,
        rng: &mut StdRng,
    ) -> SparsePoly<E::Fr, SparseTerm> {
        SparsePoly::<E::Fr, SparseTerm>::rand(degree, num_vars.unwrap(), rng)
    }

    fn rand_point<E: PairingEngine>(num_vars: Option<usize>, rng: &mut StdRng) -> Vec<E::Fr> {
        let num_vars = num_vars.unwrap();
        let mut point = Vec::with_capacity(num_vars);
        for _ in 0..num_vars {
            point.push(E::Fr::rand(rng));
        }
        point
    }

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        let num_vars = Some(10);
        single_poly_test::<_, _, PC_Bls12_377>(
            num_vars,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, PC_Bls12_381>(
            num_vars,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        let num_vars = Some(10);
        full_end_to_end_test::<_, _, PC_Bls12_377>(
            num_vars,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_test::<_, _, PC_Bls12_381>(
            num_vars,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        let num_vars = Some(10);
        single_equation_test::<_, _, PC_Bls12_377>(
            num_vars,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        single_equation_test::<_, _, PC_Bls12_381>(
            num_vars,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        let num_vars = Some(10);
        two_equation_test::<_, _, PC_Bls12_377>(
            num_vars,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_test::<_, _, PC_Bls12_381>(
            num_vars,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        let num_vars = Some(10);
        full_end_to_end_equation_test::<_, _, PC_Bls12_377>(
            num_vars,
            rand_poly::<Bls12_377>,
            rand_point::<Bls12_377>,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_equation_test::<_, _, PC_Bls12_381>(
            num_vars,
            rand_poly::<Bls12_381>,
            rand_point::<Bls12_381>,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }
}
