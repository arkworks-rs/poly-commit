use crate::Vec;
use ark_ec::group::Group;
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

pub mod data_structures;
pub use data_structures::*;

pub struct InnerProductArgPCGadget<G, C, S, SV>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>> + ToConstraintFieldGadget<ConstraintF<G>>,
    S: CryptographicSponge<ConstraintF<G>>,
    SV: CryptographicSpongeVar<ConstraintF<G>, S>,
{
    pub _affine: PhantomData<G>,
    pub _curve_var: PhantomData<C>,
    pub _sponge: PhantomData<S>,
    pub _sponge_var: PhantomData<SV>,
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
            svk_var,
            commitment_vars,
            point_var,
            value_vars,
            proof_var,
            opening_challenge_vars
        )
    )]
    pub fn succinct_check<'a>(
        cs: ConstraintSystemRef<ConstraintF<G>>,
        svk_var: &SuccinctVerifierKeyVar<G, C>,
        commitment_vars: impl IntoIterator<Item = &'a CommitmentVar<G, C>>,
        point_var: &NNFieldVar<G>,
        value_vars: impl IntoIterator<Item = &'a NNFieldVar<G>>,
        proof_var: &ProofVar<G, C>,
        opening_challenge_vars: &dyn Fn(u64) -> NNFieldVar<G>,
    ) -> Result<(Boolean<ConstraintF<G>>, SuccinctCheckPolynomialVar<G>), SynthesisError> {
        let start = cs.num_constraints();
        let d = svk_var.supported_degree;

        // `log_d` is ceil(log2 (d + 1)), which is the number of steps to compute all of the challenges
        let log_d = ark_std::log2(d + 1) as usize;

        let mut combined_commitment_var = C::zero();
        let mut combined_eval_var = NNFieldVar::<G>::zero();

        for (i, (commitment_var, value_var)) in
            commitment_vars.into_iter().zip(value_vars).enumerate()
        {
            let cur_challenge_var: NNFieldVar<G> = opening_challenge_vars((2 * i) as u64);
            // TODO: A bit hacky. May want to revert or change later on?
            // TODO: Replace combined_eval_var operations with mul_without_reduce?
            if !cur_challenge_var.is_one()?.value()? {
                combined_eval_var += &((&cur_challenge_var).mul(value_var));
                combined_commitment_var += &commitment_var
                    .comm_var
                    .scalar_mul_le((cur_challenge_var.to_bits_le()?).iter())?;
            } else {
                combined_eval_var += value_var;
                combined_commitment_var += &commitment_var.comm_var;
            }

            if let Some((degree_bound, shifted_commitment_var)) = &commitment_var.shifted_comm_var {
                let shift_var = point_var.pow_by_constant(&[(d - degree_bound) as u64])?;
                let cur_challenge_var: NNFieldVar<G> = opening_challenge_vars((2 * i + 1) as u64);
                if !cur_challenge_var.is_one()?.value()? {
                    combined_eval_var += &((&cur_challenge_var).mul((&value_var).mul(&shift_var)));
                    combined_commitment_var += &shifted_commitment_var
                        .scalar_mul_le((cur_challenge_var.to_bits_le()?).iter())?;
                } else {
                    combined_eval_var += &(&value_var).mul(&shift_var);
                    combined_commitment_var += shifted_commitment_var;
                }
            }
        }

        let combined_eval_bytes_var = combined_eval_var.to_bytes()?;

        let mut point_and_combined_eval_bytes_var = point_var.to_bytes()?;
        point_and_combined_eval_bytes_var.extend_from_slice(combined_eval_bytes_var.as_slice());

        let point_and_combined_eval_fp_vars =
            point_and_combined_eval_bytes_var.to_constraint_field()?;

        if let Some((hiding_comm_var, rand_var)) = &proof_var.hiding_var {
            let mut hiding_challenge_sponge_var =
                SV::new(ns!(cs, "hiding_challenge_sponge_var").cs());

            hiding_challenge_sponge_var
                .absorb(combined_commitment_var.to_constraint_field()?.as_slice())?;
            hiding_challenge_sponge_var
                .absorb(hiding_comm_var.to_constraint_field()?.as_slice())?;
            hiding_challenge_sponge_var.absorb(point_and_combined_eval_fp_vars.as_slice())?;

            let hiding_challenge_bits_var = hiding_challenge_sponge_var.squeeze_bits(128)?;
            combined_commitment_var += &(hiding_comm_var
                .scalar_mul_le(hiding_challenge_bits_var.iter())?
                - &(svk_var.s_var.scalar_mul_le(rand_var.iter())?));
        }

        let mut round_challenge_vars = Vec::with_capacity(log_d);

        // Challenge for each round
        let mut round_challenge_sponge_var =
            SV::new(ns!(cs, "round_challenge_sponge_var_init").cs());
        round_challenge_sponge_var
            .absorb(combined_commitment_var.to_constraint_field()?.as_slice())?;
        round_challenge_sponge_var.absorb(point_and_combined_eval_fp_vars.as_slice())?;

        // Initialize challenges
        let mut round_challenge_field_elements_and_bits = round_challenge_sponge_var
            .squeeze_nonnative_field_elements_with_sizes::<G::ScalarField>(
            &[FieldElementSize::Truncated { num_bits: 128 }],
        )?;
        let mut round_challenge_var = round_challenge_field_elements_and_bits.0.pop().unwrap();
        let mut round_challenge_bits_var = round_challenge_field_elements_and_bits.1.pop().unwrap();

        let h_prime_var = svk_var
            .h_var
            .scalar_mul_le(round_challenge_bits_var.iter())?;

        let mut round_commitment_var = combined_commitment_var
            + &h_prime_var.scalar_mul_le(combined_eval_bytes_var.to_bits_le()?.iter())?;

        for (l_var, r_var) in proof_var.l_var_vec.iter().zip(&proof_var.r_var_vec) {
            let mut round_challenge_sponge_var = SV::new(ns!(cs, "round_challenge_sponge_var").cs());

            let round_challenge_bytes_var = round_challenge_bits_var
                .chunks(8)
                .map(UInt8::<ConstraintF<G>>::from_bits_le)
                .collect::<Vec<_>>();

            round_challenge_sponge_var
                .absorb(round_challenge_bytes_var.to_constraint_field()?.as_slice());
            round_challenge_sponge_var.absorb(l_var.to_constraint_field()?.as_slice());
            round_challenge_sponge_var.absorb(r_var.to_constraint_field()?.as_slice());

            // Update challenges
            round_challenge_field_elements_and_bits = round_challenge_sponge_var
                .squeeze_nonnative_field_elements_with_sizes(&[FieldElementSize::Truncated {
                    num_bits: 128,
                }])?;
            round_challenge_var = round_challenge_field_elements_and_bits.0.pop().unwrap();
            round_challenge_bits_var = round_challenge_field_elements_and_bits.1.pop().unwrap();

            // Instead of directly computing
            // `l_var.scalar_mul_le(round_challenge_var.inverse()?.to_bits_le()?.iter())`,
            // we instead allocate a variable `temp` whose value is `l_var * round_challenge.inverse()`,
            // and then enforce that `temp * round_challenge == l_var`
            let l_times_round_ch_inv = C::new_witness(ns!(cs, "l * round_ch^{-1}"), || {
                let mut l_val = l_var.value()?;
                l_val *= round_challenge_var.value()?.inverse().unwrap();
                Ok(l_val)
            })?;
            let claimed_l_var =
                l_times_round_ch_inv.scalar_mul_le(round_challenge_bits_var.iter())?;
            claimed_l_var.enforce_equal(&l_var)?;
            round_commitment_var += &l_times_round_ch_inv;
            round_commitment_var += &(r_var.scalar_mul_le(round_challenge_bits_var.iter())?);

            round_challenge_vars.push(round_challenge_var.clone());
        }

        let check_poly_var = SuccinctCheckPolynomialVar::<G>(round_challenge_vars);
        let v_prime_var = check_poly_var.evaluate(&point_var)? * &proof_var.c_var;

        let check_commitment_elem_var = CMCommitGadget::<G, C>::commit(
            &[proof_var.final_comm_key_var.clone(), h_prime_var],
            &[proof_var.c_var.to_bits_le()?, v_prime_var.to_bits_le()?],
            None,
        )?;

        let result_var = round_commitment_var.is_eq(&check_commitment_elem_var)?;
        //println!("Succinct_check {:}", cs.num_constraints() - start);
        Ok((result_var, check_poly_var))
    }

    /*
    pub fn check<'a>(
        cs: ConstraintSystemRef<ConstraintF<G>>,
        vk_var: &VerifierKeyVar<G, C>,
        commitment_vars: impl IntoIterator<Item = &'a CommitmentVar<G, C>>,
        point_var: &NNFieldVar<G>,
        value_vars: impl IntoIterator<Item = &'a NNFieldVar<G>>,
        proof_var: &ProofVar<G, C>,
        opening_challenge_vars: &dyn Fn(u64) -> NNFieldVar<G>,
    ) -> Result<Boolean<ConstraintF<G>>, SynthesisError> {
        let mut check_result_var = Boolean::TRUE;

        let d = vk_var.supported_degree();

        // `log_d` is ceil(log2 (d + 1)), which is the number of steps to compute all of the challenges
        let log_d = ark_std::log2(d + 1) as usize;

        if proof_var.l_var_vec.len() != proof_var.r_var_vec.len()
            || proof_var.l_var_vec.len() != log_d
        {
            return Ok(Boolean::FALSE);
        }

        let (succinct_check_result_var, check_poly_var) = Self::succinct_check(
            ns!(cs, "succinct_check").cs(),
            vk_var,
            commitment_vars,
            point_var,
            value_vars,
            proof_var,
            opening_challenge_vars,
        )?;

        check_result_var = check_result_var.and(&succinct_check_result_var)?;

        let check_poly_coeffs = check_poly_var.compute_coeff_vars();
        let final_key_var = CMCommitGadget::<G, C>::commit(
            vk_var.comm_key_var.as_slice(),
            check_poly_coeffs
                .iter()
                .map(|c| c.to_bits_le())
                .collect::<Result<Vec<_>, SynthesisError>>()?
                .as_slice(),
            None,
        )?;

        check_result_var =
            check_result_var.and(&(final_key_var.is_eq(&proof_var.final_comm_key_var)?))?;
        Ok(check_result_var)
    }

     */
}

#[cfg(test)]
pub mod tests {
    use crate::ipa_pc::constraints::{
        CommitmentVar, InnerProductArgPCGadget, NNFieldVar, ProofVar, SuccinctVerifierKeyVar,
        VerifierKeyVar,
    };
    use crate::ipa_pc::{InnerProductArgPC, SuccinctVerifierKey};
    use crate::{LabeledPolynomial, PolynomialCommitment, PolynomialLabel};
    use ark_ed_on_bls12_381::constraints::EdwardsVar;
    use ark_ed_on_bls12_381::EdwardsAffine;
    use ark_ff::{test_rng, One, UniformRand};
    use ark_poly::polynomial::univariate::DensePolynomial;
    use ark_poly::{univariate::DensePolynomial as DensePoly, Polynomial, UVPolynomial};
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
        let vk_var: SuccinctVerifierKeyVar<G, C> =
            SuccinctVerifierKeyVar::<G, C>::new_constant(cs.clone(), svk).unwrap();

        let commitment_var =
            CommitmentVar::<G, C>::new_constant(cs.clone(), commitment[0].clone()).unwrap();
        let point_var = NNFieldVar::<G>::new_constant(cs.clone(), point).unwrap();
        let value_var = NNFieldVar::<G>::new_constant(cs.clone(), value).unwrap();
        let proof_var = ProofVar::<G, C>::new_constant(cs.clone(), proof).unwrap();

        /*
        let check = InnerProductArgPCGadget::<G, C>::succinct_check(
            cs.clone(),
            &vk_var,
            vec![&commitment_var],
            &point_var,
            vec![&value_var],
            &proof_var,
            &|_| NNFieldVar::<G>::one(),
        )
        .unwrap();
        check.0.enforce_equal(&Boolean::TRUE);

         */

        let check = InnerProductArgPCGadget::<G, C, PoseidonSpongeVar<CF>>::succinct_check(
            cs.clone(),
            &vk_var,
            vec![&commitment_var],
            &point_var,
            vec![&value_var],
            &proof_var,
            &|_| NNFieldVar::<G>::one(),
        )
        .unwrap();

        check.0.enforce_equal(&Boolean::TRUE);

        assert!(cs.is_satisfied().unwrap());
    }
}
