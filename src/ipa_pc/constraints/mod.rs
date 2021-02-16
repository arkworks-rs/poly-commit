use crate::Vec;
use ark_ec::AffineCurve;
use ark_ff::Field;
use ark_r1cs_std::bits::boolean::Boolean;
use ark_r1cs_std::bits::uint8::UInt8;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::groups::CurveVar;
use ark_r1cs_std::{R1CSVar, ToConstraintFieldGadget};
use ark_r1cs_std::{ToBitsGadget, ToBytesGadget};
use ark_relations::ns;
use ark_relations::r1cs::{ConstraintSystemRef, SynthesisError};
use ark_sponge::constraints::CryptographicSpongeVar;
use ark_sponge::{FieldElementSize, CryptographicSponge};
use ark_std::marker::PhantomData;
use ark_std::ops::Mul;

/// Gadget versions of data structures for PC_IPA
pub mod data_structures;
pub use data_structures::*;

/// Wraps some common gadgets for PC+IPA.
pub struct InnerProductArgPCGadget<G, C, S, SV>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>> + ToConstraintFieldGadget<ConstraintF<G>>,
    S: CryptographicSponge<ConstraintF<G>>,
    SV: CryptographicSpongeVar<ConstraintF<G>, S>,
{
    _affine: PhantomData<G>,
    _curve_var: PhantomData<C>,
    _sponge: PhantomData<S>,
    _sponge_var: PhantomData<SV>,
}

impl<G, C, S, SV> InnerProductArgPCGadget<G, C, S, SV>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>> + ToConstraintFieldGadget<ConstraintF<G>>,
    S: CryptographicSponge<ConstraintF<G>>,
    SV: CryptographicSpongeVar<ConstraintF<G>, S>,
{
    /// The succinct portion of `PC::check`. This algorithm runs in time
    /// O(log d), where d is the degree of the committed polynomials.
    #[tracing::instrument(
        target = "r1cs",
        skip(
            cs,
            svk,
            commitments,
            point,
            values,
            proof,
            opening_challenges
        )
    )]
    pub fn succinct_check<'a>(
        cs: ConstraintSystemRef<ConstraintF<G>>,
        svk: &SuccinctVerifierKeyVar<G, C>,
        commitments: impl IntoIterator<Item = &'a CommitmentVar<G, C>>,
        point: &NNFieldVar<G>,
        values: impl IntoIterator<Item = &'a NNFieldVar<G>>,
        proof: &ProofVar<G, C>,
        opening_challenges: &dyn Fn(u64) -> NNFieldVar<G>,
    ) -> Result<(Boolean<ConstraintF<G>>, SuccinctCheckPolynomialVar<G>), SynthesisError> {
        let d = svk.supported_degree;

        // `log_d` is ceil(log2 (d + 1)), which is the number of steps to compute all of the challenges
        let log_d = ark_std::log2(d + 1) as usize;

        let mut combined_commitment = C::zero();
        let mut combined_eval = NNFieldVar::<G>::zero();

        for (i, (commitment, value)) in
            commitments.into_iter().zip(values).enumerate()
        {
            let cur_challenge: NNFieldVar<G> = opening_challenges((2 * i) as u64);
            // TODO: A bit hacky. May want to revert or change later on?
            // TODO: Replace combined_eval operations with mul_without_reduce?
            if !cur_challenge.is_one()?.value()? {
                combined_eval += &((&cur_challenge).mul(value));
                combined_commitment += &commitment
                    .comm
                    .scalar_mul_le((cur_challenge.to_bits_le()?).iter())?;
            } else {
                combined_eval += value;
                combined_commitment += &commitment.comm;
            }

            if let Some((degree_bound, shifted_commitment)) = &commitment.shifted_comm {
                let shift = point.pow_by_constant(&[(d - degree_bound) as u64])?;
                let cur_challenge: NNFieldVar<G> = opening_challenges((2 * i + 1) as u64);
                if !cur_challenge.is_one()?.value()? {
                    combined_eval += &((&cur_challenge).mul((&value).mul(&shift)));
                    combined_commitment += &shifted_commitment
                        .scalar_mul_le((cur_challenge.to_bits_le()?).iter())?;
                } else {
                    combined_eval += &(&value).mul(&shift);
                    combined_commitment += shifted_commitment;
                }
            }
        }

        let combined_eval_bytes = combined_eval.to_bytes()?;

        let mut point_and_combined_eval_bytes = point.to_bytes()?;
        point_and_combined_eval_bytes.extend_from_slice(combined_eval_bytes.as_slice());

        let point_and_combined_eval_fps =
            point_and_combined_eval_bytes.to_constraint_field()?;

        if let Some((hiding_comm, rand)) = &proof.hiding {
            let mut hiding_challenge_sponge =
                SV::new(ns!(cs, "hiding_challenge_sponge").cs());

            hiding_challenge_sponge
                .absorb(combined_commitment.to_constraint_field()?.as_slice())?;
            hiding_challenge_sponge
                .absorb(hiding_comm.to_constraint_field()?.as_slice())?;
            hiding_challenge_sponge.absorb(point_and_combined_eval_fps.as_slice())?;

            let hiding_challenge_bits = hiding_challenge_sponge.squeeze_bits(128)?;
            combined_commitment += &(hiding_comm
                .scalar_mul_le(hiding_challenge_bits.iter())?
                - &(svk.s.scalar_mul_le(rand.iter())?));
        }

        let mut round_challenges = Vec::with_capacity(log_d);

        // Challenge for each round
        let mut round_challenge_sponge =
            SV::new(ns!(cs, "round_challenge_sponge_init").cs());
        round_challenge_sponge
            .absorb(combined_commitment.to_constraint_field()?.as_slice())?;
        round_challenge_sponge.absorb(point_and_combined_eval_fps.as_slice())?;

        // Initialize challenges
        let mut round_challenge_field_elements_and_bits = round_challenge_sponge
            .squeeze_nonnative_field_elements_with_sizes::<G::ScalarField>(
            &[FieldElementSize::Truncated { num_bits: 128 }],
        )?;
        let mut round_challenge;
        let mut round_challenge_bits = round_challenge_field_elements_and_bits.1.pop().unwrap();

        let h_prime = svk
            .h
            .scalar_mul_le(round_challenge_bits.iter())?;

        let mut round_commitment = combined_commitment
            + &h_prime.scalar_mul_le(combined_eval_bytes.to_bits_le()?.iter())?;

        for (l, r) in proof.l_vec.iter().zip(&proof.r_vec) {
            let mut round_challenge_sponge = SV::new(ns!(cs, "round_challenge_sponge").cs());

            let round_challenge_bytes = round_challenge_bits
                .chunks(8)
                .map(UInt8::<ConstraintF<G>>::from_bits_le)
                .collect::<Vec<_>>();

            round_challenge_sponge
                .absorb(round_challenge_bytes.to_constraint_field()?.as_slice())?;
            round_challenge_sponge.absorb(l.to_constraint_field()?.as_slice())?;
            round_challenge_sponge.absorb(r.to_constraint_field()?.as_slice())?;

            // Update challenges
            round_challenge_field_elements_and_bits = round_challenge_sponge
                .squeeze_nonnative_field_elements_with_sizes(&[FieldElementSize::Truncated {
                    num_bits: 128,
                }])?;
            round_challenge = round_challenge_field_elements_and_bits.0.pop().unwrap();
            round_challenge_bits = round_challenge_field_elements_and_bits.1.pop().unwrap();

            // Instead of directly computing
            // `l.scalar_mul_le(round_challenge.inverse()?.to_bits_le()?.iter())`,
            // we instead allocate a variable `temp` whose value is `l * round_challenge.inverse()`,
            // and then enforce that `temp * round_challenge == l`
            let l_times_round_ch_inv = C::new_witness(ns!(cs, "l * round_ch^{-1}"), || {
                let mut l_val = l.value()?;
                l_val *= round_challenge.value()?.inverse().unwrap();
                Ok(l_val)
            })?;
            let claimed_l =
                l_times_round_ch_inv.scalar_mul_le(round_challenge_bits.iter())?;
            claimed_l.enforce_equal(&l)?;
            round_commitment += &l_times_round_ch_inv;
            round_commitment += &(r.scalar_mul_le(round_challenge_bits.iter())?);

            round_challenges.push(round_challenge.clone());
        }

        let check_poly = SuccinctCheckPolynomialVar::<G>(round_challenges);
        let v_prime = check_poly.evaluate(&point)? * &proof.c;

        let check_commitment_elem = CMCommitGadget::<G, C>::commit(
            &[proof.final_comm_key.clone(), h_prime],
            &[proof.c.to_bits_le()?, v_prime.to_bits_le()?],
            None,
        )?;

        let result = round_commitment.is_eq(&check_commitment_elem)?;
        //println!("Succinct_check {:}", cs.num_constraints() - start);
        Ok((result, check_poly))
    }

    /*
    pub fn check<'a>(
        cs: ConstraintSystemRef<ConstraintF<G>>,
        vk: &VerifierKeyVar<G, C>,
        commitments: impl IntoIterator<Item = &'a CommitmentVar<G, C>>,
        point: &NNFieldVar<G>,
        values: impl IntoIterator<Item = &'a NNFieldVar<G>>,
        proof: &ProofVar<G, C>,
        opening_challenges: &dyn Fn(u64) -> NNFieldVar<G>,
    ) -> Result<Boolean<ConstraintF<G>>, SynthesisError> {
        let mut check_result = Boolean::TRUE;

        let d = vk.supported_degree();

        // `log_d` is ceil(log2 (d + 1)), which is the number of steps to compute all of the challenges
        let log_d = ark_std::log2(d + 1) as usize;

        if proof.l_vec.len() != proof.r_vec.len()
            || proof.l_vec.len() != log_d
        {
            return Ok(Boolean::FALSE);
        }

        let (succinct_check_result, check_poly) = Self::succinct_check(
            ns!(cs, "succinct_check").cs(),
            vk,
            commitments,
            point,
            values,
            proof,
            opening_challenges,
        )?;

        check_result = check_result.and(&succinct_check_result)?;

        let check_poly_coeffs = check_poly.compute_coeffs();
        let final_key = CMCommitGadget::<G, C>::commit(
            vk.comm_key.as_slice(),
            check_poly_coeffs
                .iter()
                .map(|c| c.to_bits_le())
                .collect::<Result<Vec<_>, SynthesisError>>()?
                .as_slice(),
            None,
        )?;

        check_result =
            check_result.and(&(final_key.is_eq(&proof.final_comm_key)?))?;
        Ok(check_result)
    }

     */
}

#[cfg(test)]
pub mod tests {
    use crate::ipa_pc::constraints::{
        CommitmentVar, InnerProductArgPCGadget, NNFieldVar, ProofVar, SuccinctVerifierKeyVar,
    };
    use crate::ipa_pc::{InnerProductArgPC, SuccinctVerifierKey};
    use crate::{LabeledPolynomial, PolynomialCommitment, PolynomialLabel};
    use ark_std::{test_rng, UniformRand};
    use ark_ff::One;
    use ark_poly::polynomial::univariate::DensePolynomial;
    use ark_poly::{univariate::DensePolynomial as DensePoly, UVPolynomial};
    use ark_r1cs_std::alloc::AllocVar;
    use ark_r1cs_std::bits::boolean::Boolean;
    use ark_r1cs_std::eq::EqGadget;
    use ark_r1cs_std::fields::FieldVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_sponge::poseidon::constraints::PoseidonSpongeVar;
    use ark_sponge::poseidon::PoseidonSponge;
    use blake2::Blake2s;

    type G = ark_pallas::Affine;
    type C = ark_pallas::constraints::GVar;
    type F = ark_pallas::Fr;
    type CF = ark_pallas::Fq;

    type UniPoly = DensePoly<F>;
    type PC<E, D, P, CF, S> = InnerProductArgPC<E, D, P, CF, S>;
    type PC_JJB2S = PC<G, Blake2s, UniPoly, CF, PoseidonSponge<CF>>;

    #[test]
    pub fn basic() {
        let mut rng = test_rng();
        let random_polynomial = DensePolynomial::rand(20, &mut rng);
        let labeled_random_polynomial = LabeledPolynomial::new(
            PolynomialLabel::new(),
            random_polynomial,
            Some(20),
            Some(20),
        );

        let pp = PC_JJB2S::setup(20, None, &mut rng).unwrap();
        let (ck, vk) = PC_JJB2S::trim(&pp, 20, 0, None).unwrap();
        let (commitment, randomness) =
            PC_JJB2S::commit(&ck, vec![&labeled_random_polynomial], Some(&mut rng)).unwrap();

        let point = F::rand(&mut rng);
        let value = labeled_random_polynomial.evaluate(&point);
        let proof = PC_JJB2S::open(
            &ck,
            vec![&labeled_random_polynomial],
            &commitment,
            &point,
            F::one(),
            &randomness,
            Some(&mut rng),
        )
        .unwrap();

        assert!(PC_JJB2S::check(
            &vk,
            &commitment,
            &point,
            vec![value],
            &proof,
            F::one(),
            Some(&mut rng)
        )
        .unwrap());

        let cs = ConstraintSystem::<CF>::new_ref();
        let svk = SuccinctVerifierKey::from_vk(&vk);
        let vk: SuccinctVerifierKeyVar<G, C> =
            SuccinctVerifierKeyVar::<G, C>::new_constant(cs.clone(), svk).unwrap();

        let commitment =
            CommitmentVar::<G, C>::new_constant(cs.clone(), commitment[0].clone()).unwrap();
        let point = NNFieldVar::<G>::new_constant(cs.clone(), point).unwrap();
        let value = NNFieldVar::<G>::new_constant(cs.clone(), value).unwrap();
        let proof = ProofVar::<G, C>::new_constant(cs.clone(), proof).unwrap();

        /*
        let check = InnerProductArgPCGadget::<G, C>::succinct_check(
            cs.clone(),
            &vk,
            vec![&commitment],
            &point,
            vec![&value],
            &proof,
            &|_| NNFieldVar::<G>::one(),
        )
        .unwrap();
        check.0.enforce_equal(&Boolean::TRUE);

         */

        let check = InnerProductArgPCGadget::<G, C, PoseidonSpongeVar<CF>>::succinct_check(
            cs.clone(),
            &vk,
            vec![&commitment],
            &point,
            vec![&value],
            &proof,
            &|_| NNFieldVar::<G>::one(),
        )
        .unwrap();

        check.0.enforce_equal(&Boolean::TRUE).unwrap();

        assert!(cs.is_satisfied().unwrap());
    }
}
