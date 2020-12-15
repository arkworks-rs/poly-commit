use crate::Vec;
use ark_ec::AffineCurve;
use ark_marlin::fiat_shamir::constraints::FiatShamirRngVar;
use ark_marlin::fiat_shamir::FiatShamirRng;
use ark_r1cs_std::bits::boolean::Boolean;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::groups::CurveVar;
use ark_r1cs_std::{ToBitsGadget, ToBytesGadget};
use ark_relations::r1cs::{ConstraintSystemRef, SynthesisError};
use ark_std::marker::PhantomData;
use ark_std::ops::Mul;

pub mod data_structures;
pub use data_structures::*;

pub struct InnerProductArgPCGadget<G, C, S, SV>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
    S: FiatShamirRng<G::ScalarField, ConstraintF<G>>,
    SV: FiatShamirRngVar<G::ScalarField, ConstraintF<G>, S>,
{
    pub _affine: PhantomData<G>,
    pub _curve_var: PhantomData<C>,
    pub _sponge: PhantomData<S>,
    pub _sponge_var: PhantomData<SV>,
}

impl<G, C, S, SV> InnerProductArgPCGadget<G, C, S, SV>
where
    G: AffineCurve,
    C: CurveVar<G::Projective, ConstraintF<G>>,
    S: FiatShamirRng<G::ScalarField, ConstraintF<G>>,
    SV: FiatShamirRngVar<G::ScalarField, ConstraintF<G>, S>,
{
    /// The succinct portion of `PC::check`. This algorithm runs in time
    /// O(log d), where d is the degree of the committed polynomials.
    pub fn succinct_check<'a>(
        cs: ConstraintSystemRef<ConstraintF<G>>,
        vk_var: &VerifierKeyVar<G, C>,
        commitment_vars: impl IntoIterator<Item = &'a CommitmentVar<G, C>>,
        point_var: &NNFieldVar<G>,
        value_vars: impl IntoIterator<Item = &'a NNFieldVar<G>>,
        proof_var: &ProofVar<G, C>,
        opening_challenge_vars: &dyn Fn(u64) -> NNFieldVar<G>,
    ) -> Result<(Boolean<ConstraintF<G>>, SuccinctCheckPolynomialVar<G>), SynthesisError> {
        // To substitute sponge
        // TODO: Remove when implemented sponge
        //let two = G::ScalarField::one() + &G::ScalarField::one();
        //let two_var = NNFieldVar::<G>::new_constant(cs.clone(), two)?;

        let d = vk_var.supported_degree();

        // `log_d` is ceil(log2 (d + 1)), which is the number of steps to compute all of the challenges
        let log_d = ark_std::log2(d + 1) as usize;

        let mut combined_commitment_var = C::zero();
        let mut combined_eval_var = NNFieldVar::<G>::zero();

        for (i, (commitment_var, value_var)) in
            commitment_vars.into_iter().zip(value_vars).enumerate()
        {
            let cur_challenge_var: NNFieldVar<G> = opening_challenge_vars((2 * i) as u64);
            combined_eval_var += &((&cur_challenge_var).mul(value_var));
            combined_commitment_var += &commitment_var
                .comm_var
                .scalar_mul_le((cur_challenge_var.to_bits_le()?).iter())?;

            if let Some((degree_bound, shifted_commitment_var)) = &commitment_var.shifted_comm_var {
                let cur_challenge_var: NNFieldVar<G> = opening_challenge_vars((2 * i + 1) as u64);
                let shift_var = point_var.pow_by_constant(&[(d - degree_bound) as u64])?;
                combined_eval_var += &((&cur_challenge_var).mul((&value_var).mul(&shift_var)));
                combined_commitment_var += &shifted_commitment_var
                    .scalar_mul_le((cur_challenge_var.to_bits_le()?).iter())?;
            }
        }

        if let Some((hiding_comm_var, rand_var)) = &proof_var.hiding_var {
            let mut hiding_challenge_sponge_var = SV::new(cs.clone());
            hiding_challenge_sponge_var
                .absorb_bytes(combined_commitment_var.to_bytes()?.as_slice())?;
            hiding_challenge_sponge_var.absorb_bytes(hiding_comm_var.to_bytes()?.as_slice())?;
            hiding_challenge_sponge_var.absorb_bytes(point_var.to_bytes()?.as_slice())?;
            hiding_challenge_sponge_var.absorb_bytes(combined_eval_var.to_bytes()?.as_slice())?;

            let hiding_challenge_var = hiding_challenge_sponge_var
                .squeeze_field_elements(1)?
                .pop()
                .unwrap();

            combined_commitment_var += &(hiding_comm_var
                .scalar_mul_le(hiding_challenge_var.to_bits_le()?.iter())?
                - &(vk_var.s_var.scalar_mul_le(rand_var.to_bits_le()?.iter())?));
        }

        let mut round_challenge_vars = Vec::with_capacity(log_d);

        // Challenge for each round
        let mut round_challenge_sponge_var = SV::new(cs.clone());
        round_challenge_sponge_var.absorb_bytes(combined_commitment_var.to_bytes()?.as_slice())?;
        round_challenge_sponge_var.absorb_bytes(point_var.to_bytes()?.as_slice());
        round_challenge_sponge_var.absorb_bytes(combined_eval_var.to_bytes()?.as_slice());

        let mut round_challenge_var = round_challenge_sponge_var
            .squeeze_field_elements(1)?
            .pop()
            .unwrap();

        let h_prime_var = vk_var
            .h_var
            .scalar_mul_le(round_challenge_var.to_bits_le()?.iter())?;

        let mut round_commitment_var = combined_commitment_var
            + &h_prime_var.scalar_mul_le(combined_eval_var.to_bits_le()?.iter())?;

        for (l_var, r_var) in proof_var.l_var_vec.iter().zip(&proof_var.r_var_vec) {
            let mut round_challenge_sponge_var = SV::new(cs.clone());
            round_challenge_sponge_var.absorb_bytes(round_challenge_var.to_bytes()?.as_slice())?;
            round_challenge_sponge_var.absorb_bytes(l_var.to_bytes()?.as_slice())?;
            round_challenge_sponge_var.absorb_bytes(r_var.to_bytes()?.as_slice())?;

            round_challenge_var = round_challenge_sponge_var
                .squeeze_field_elements(1)?
                .pop()
                .unwrap();

            round_commitment_var +=
                &(l_var.scalar_mul_le(round_challenge_var.inverse()?.to_bits_le()?.iter()))?;
            round_commitment_var +=
                &(r_var.scalar_mul_le(round_challenge_var.to_bits_le()?.iter())?);

            round_challenge_vars.push(round_challenge_var.clone());
        }

        let check_poly_var = SuccinctCheckPolynomialVar::<G>(round_challenge_vars);
        let v_prime_var = check_poly_var.evaluate(&point_var)? * &proof_var.c_var;

        let check_commitment_elem_var = CMCommitGadget::<G, C>::commit(
            &[proof_var.final_comm_key_var.clone(), h_prime_var],
            &[proof_var.c_var.clone(), v_prime_var],
            None,
        )?;

        let result_var = round_commitment_var.is_eq(&check_commitment_elem_var)?;
        Ok((result_var, check_poly_var))
    }

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
            cs.clone(),
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
            check_poly_coeffs.as_slice(),
            None,
        )?;

        check_result_var =
            check_result_var.and(&(final_key_var.is_eq(&proof_var.final_comm_key_var)?))?;
        Ok(check_result_var)
    }
}

#[cfg(test)]
pub mod tests {
    use crate::ipa_pc::constraints::{
        CommitmentVar, InnerProductArgPCGadget, NNFieldVar, ProofVar, VerifierKeyVar,
    };
    use crate::ipa_pc::InnerProductArgPC;
    use crate::{LabeledPolynomial, PolynomialCommitment, PolynomialLabel};
    use ark_ed_on_bls12_381::constraints::EdwardsVar;
    use ark_ed_on_bls12_381::EdwardsAffine;
    use ark_ff::{test_rng, One, UniformRand};
    use ark_marlin::fiat_shamir::constraints::FiatShamirAlgebraicSpongeRngVar;
    use ark_marlin::fiat_shamir::poseidon::constraints::PoseidonSpongeVar;
    use ark_marlin::fiat_shamir::poseidon::PoseidonSponge;
    use ark_marlin::fiat_shamir::FiatShamirAlgebraicSpongeRng;
    use ark_poly::polynomial::univariate::DensePolynomial;
    use ark_poly::{univariate::DensePolynomial as DensePoly, Polynomial, UVPolynomial};
    use ark_r1cs_std::alloc::AllocVar;
    use ark_r1cs_std::bits::boolean::Boolean;
    use ark_r1cs_std::eq::EqGadget;
    use ark_r1cs_std::fields::FieldVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_sponge::poseidon::PoseidonSpongeWrapper;
    use blake2::Blake2s;

    type G = EdwardsAffine;
    type C = EdwardsVar;
    type F = ark_ed_on_bls12_381::Fr;
    type ConstraintF = ark_ed_on_bls12_381::Fq;

    type UniPoly = DensePoly<F>;
    type PC<E, D, P, S> = InnerProductArgPC<E, D, P, S>;
    type PC_JJB2S = PC<G, Blake2s, UniPoly, PoseidonSpongeWrapper<F, ConstraintF>>;

    type Poseidon = FiatShamirAlgebraicSpongeRng<F, ConstraintF, PoseidonSponge<ConstraintF>>;
    type PoseidonVar = FiatShamirAlgebraicSpongeRngVar<
        F,
        ConstraintF,
        PoseidonSponge<ConstraintF>,
        PoseidonSpongeVar<ConstraintF>,
    >;

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

        let cs = ConstraintSystem::<ConstraintF>::new_ref();
        let vk_var: VerifierKeyVar<G, C> =
            VerifierKeyVar::<G, C>::new_constant(cs.clone(), vk).unwrap();
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

        let check = InnerProductArgPCGadget::<G, C, Poseidon, PoseidonVar>::check(
            cs.clone(),
            &vk_var,
            vec![&commitment_var],
            &point_var,
            vec![&value_var],
            &proof_var,
            &|_| NNFieldVar::<G>::one(),
        )
        .unwrap();

        check.enforce_equal(&Boolean::TRUE);

        assert!(cs.is_satisfied().unwrap());
    }
}
