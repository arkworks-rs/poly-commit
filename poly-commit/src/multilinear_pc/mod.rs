use crate::multilinear_pc::data_structures::{
    Commitment, CommitterKey, Proof, UniversalParams, VerifierKey,
};
use ark_ec::AffineRepr;
use ark_ec::{pairing::Pairing, CurveGroup};
use ark_ec::{scalar_mul::fixed_base::FixedBase, VariableBaseMSM};
use ark_ff::{Field, PrimeField};
use ark_ff::{One, Zero};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::collections::LinkedList;
use ark_std::iter::FromIterator;
use ark_std::marker::PhantomData;
use ark_std::ops::Mul;
use ark_std::rand::RngCore;
use ark_std::vec::Vec;
use ark_std::UniformRand;

/// data structures used by multilinear extension commitment scheme
pub mod data_structures;

/// Polynomial Commitment Scheme on multilinear extensions.
pub struct MultilinearPC<E: Pairing> {
    _engine: PhantomData<E>,
}

impl<E: Pairing> MultilinearPC<E> {
    /// setup
    pub fn setup<R: RngCore>(num_vars: usize, rng: &mut R) -> UniversalParams<E> {
        assert!(num_vars > 0, "constant polynomial not supported");
        let g: E::G1 = E::G1::rand(rng);
        let h: E::G2 = E::G2::rand(rng);
        let g = g.into_affine();
        let h = h.into_affine();
        let mut powers_of_g = Vec::new();
        let mut powers_of_h = Vec::new();
        let t: Vec<_> = (0..num_vars).map(|_| E::ScalarField::rand(rng)).collect();
        let scalar_bits = E::ScalarField::MODULUS_BIT_SIZE as usize;

        let mut eq: LinkedList<DenseMultilinearExtension<E::ScalarField>> =
            LinkedList::from_iter(eq_extension(&t).into_iter());
        let mut eq_arr = LinkedList::new();
        let mut base = eq.pop_back().unwrap().evaluations;

        for i in (0..num_vars).rev() {
            eq_arr.push_front(remove_dummy_variable(&base, i));
            if i != 0 {
                let mul = eq.pop_back().unwrap().evaluations;
                base = base
                    .into_iter()
                    .zip(mul.into_iter())
                    .map(|(a, b)| a * &b)
                    .collect();
            }
        }

        let mut pp_powers = Vec::new();
        let mut total_scalars = 0;
        for i in 0..num_vars {
            let eq = eq_arr.pop_front().unwrap();
            let pp_k_powers = (0..(1 << (num_vars - i))).map(|x| eq[x]);
            pp_powers.extend(pp_k_powers);
            total_scalars += 1 << (num_vars - i);
        }
        let window_size = FixedBase::get_mul_window_size(total_scalars);
        let g_table = FixedBase::get_window_table(scalar_bits, window_size, g.into_group());
        let h_table = FixedBase::get_window_table(scalar_bits, window_size, h.into_group());

        let pp_g = E::G1::normalize_batch(&FixedBase::msm(
            scalar_bits,
            window_size,
            &g_table,
            &pp_powers,
        ));
        let pp_h = E::G2::normalize_batch(&FixedBase::msm(
            scalar_bits,
            window_size,
            &h_table,
            &pp_powers,
        ));
        let mut start = 0;
        for i in 0..num_vars {
            let size = 1 << (num_vars - i);
            let pp_k_g = (&pp_g[start..(start + size)]).to_vec();
            let pp_k_h = (&pp_h[start..(start + size)]).to_vec();
            powers_of_g.push(pp_k_g);
            powers_of_h.push(pp_k_h);
            start += size;
        }

        // uncomment to measure the time for calculating vp
        // let vp_generation_timer = start_timer!(|| "VP generation");
        let g_mask = {
            let window_size = FixedBase::get_mul_window_size(num_vars);
            let g_table = FixedBase::get_window_table(scalar_bits, window_size, g.into_group());
            E::G1::normalize_batch(&FixedBase::msm(scalar_bits, window_size, &g_table, &t))
        };
        // end_timer!(vp_generation_timer);

        UniversalParams {
            num_vars,
            g,
            g_mask,
            h,
            powers_of_g,
            powers_of_h,
        }
    }

    /// Trim the universal parameters to specialize the public parameters
    /// for multilinear polynomials to the given `supported_num_vars`, and returns committer key and verifier key.
    /// `supported_num_vars` should be in range `1..=params.num_vars`
    pub fn trim(
        params: &UniversalParams<E>,
        supported_num_vars: usize,
    ) -> (CommitterKey<E>, VerifierKey<E>) {
        assert!(supported_num_vars <= params.num_vars);
        let to_reduce = params.num_vars - supported_num_vars;
        let ck = CommitterKey {
            powers_of_h: (&params.powers_of_h[to_reduce..]).to_vec(),
            powers_of_g: (&params.powers_of_g[to_reduce..]).to_vec(),
            g: params.g,
            h: params.h,
            nv: supported_num_vars,
        };
        let vk = VerifierKey {
            nv: supported_num_vars,
            g: params.g,
            h: params.h,
            g_mask_random: (&params.g_mask[to_reduce..]).to_vec(),
        };
        (ck, vk)
    }

    /// commit
    pub fn commit(
        ck: &CommitterKey<E>,
        polynomial: &impl MultilinearExtension<E::ScalarField>,
    ) -> Commitment<E> {
        let nv = polynomial.num_vars();
        let scalars: Vec<_> = polynomial
            .to_evaluations()
            .into_iter()
            .map(|x| x.into_bigint())
            .collect();
        let g_product =
            <E::G1 as VariableBaseMSM>::msm_bigint(&ck.powers_of_g[0], scalars.as_slice())
                .into_affine();
        Commitment { nv, g_product }
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the same.
    pub fn open(
        ck: &CommitterKey<E>,
        polynomial: &impl MultilinearExtension<E::ScalarField>,
        point: &[E::ScalarField],
    ) -> Proof<E> {
        assert_eq!(polynomial.num_vars(), ck.nv, "Invalid size of polynomial");
        let nv = polynomial.num_vars();
        let mut r: Vec<Vec<E::ScalarField>> = (0..nv + 1).map(|_| Vec::new()).collect();
        let mut q: Vec<Vec<E::ScalarField>> = (0..nv + 1).map(|_| Vec::new()).collect();

        r[nv] = polynomial.to_evaluations();

        let mut proofs = Vec::new();
        for i in 0..nv {
            let k = nv - i;
            let point_at_k = point[i];
            q[k] = (0..(1 << (k - 1)))
                .map(|_| E::ScalarField::zero())
                .collect();
            r[k - 1] = (0..(1 << (k - 1)))
                .map(|_| E::ScalarField::zero())
                .collect();
            for b in 0..(1 << (k - 1)) {
                q[k][b] = r[k][(b << 1) + 1] - &r[k][b << 1];
                r[k - 1][b] = r[k][b << 1] * &(E::ScalarField::one() - &point_at_k)
                    + &(r[k][(b << 1) + 1] * &point_at_k);
            }
            let scalars: Vec<_> = (0..(1 << k))
                .map(|x| q[k][x >> 1].into_bigint()) // fine
                .collect();

            let pi_h =
                <E::G2 as VariableBaseMSM>::msm_bigint(&ck.powers_of_h[i], &scalars).into_affine(); // no need to move outside and partition
            proofs.push(pi_h);
        }

        Proof { proofs }
    }

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    pub fn check<'a>(
        vk: &VerifierKey<E>,
        commitment: &Commitment<E>,
        point: &[E::ScalarField],
        value: E::ScalarField,
        proof: &Proof<E>,
    ) -> bool {
        let left = E::pairing(commitment.g_product.into_group() - &vk.g.mul(value), vk.h);

        let scalar_size = E::ScalarField::MODULUS_BIT_SIZE as usize;
        let window_size = FixedBase::get_mul_window_size(vk.nv);

        let g_table = FixedBase::get_window_table(scalar_size, window_size, vk.g.into_group());
        let g_mul: Vec<E::G1> = FixedBase::msm(scalar_size, window_size, &g_table, point);

        let pairing_lefts: Vec<_> = (0..vk.nv)
            .map(|i| vk.g_mask_random[i].into_group() - &g_mul[i])
            .collect();
        let pairing_lefts: Vec<E::G1Affine> = E::G1::normalize_batch(&pairing_lefts);
        let pairing_lefts: Vec<E::G1Prepared> = pairing_lefts
            .into_iter()
            .map(|x| E::G1Prepared::from(x))
            .collect();

        let pairing_rights: Vec<E::G2Prepared> = proof
            .proofs
            .iter()
            .map(|x| E::G2Prepared::from(*x))
            .collect();

        let right = E::multi_pairing(pairing_lefts, pairing_rights);
        left == right
    }
}

/// fix first `pad` variables of `poly` represented in evaluation form to zero
fn remove_dummy_variable<F: Field>(poly: &[F], pad: usize) -> Vec<F> {
    if pad == 0 {
        return poly.to_vec();
    }
    if !poly.len().is_power_of_two() {
        panic!("Size of polynomial should be power of two. ")
    }
    let nv = ark_std::log2(poly.len()) as usize - pad;
    let table: Vec<_> = (0..(1 << nv)).map(|x| poly[x << pad]).collect();
    table
}

/// generate eq(t,x), a product of multilinear polynomials with fixed t.
/// eq(a,b) is takes extensions of a,b in {0,1}^num_vars such that if a and b in {0,1}^num_vars are equal
/// then this polynomial evaluates to 1.
fn eq_extension<F: Field>(t: &[F]) -> Vec<DenseMultilinearExtension<F>> {
    let dim = t.len();
    let mut result = Vec::new();
    for i in 0..dim {
        let mut poly = Vec::with_capacity(1 << dim);
        for x in 0..(1 << dim) {
            let xi = if x >> i & 1 == 1 { F::one() } else { F::zero() };
            let ti = t[i];
            let ti_xi = ti * xi;
            poly.push(ti_xi + ti_xi - xi - ti + F::one());
        }
        result.push(DenseMultilinearExtension::from_evaluations_vec(dim, poly));
    }

    result
}

#[cfg(test)]
mod tests {
    use crate::ark_std::UniformRand;
    use crate::multilinear_pc::data_structures::UniversalParams;
    use crate::multilinear_pc::MultilinearPC;
    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::Pairing;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};
    use ark_std::rand::RngCore;
    use ark_std::test_rng;
    use ark_std::vec::Vec;
    type E = Bls12_381;
    type Fr = <E as Pairing>::ScalarField;

    fn test_polynomial<R: RngCore>(
        uni_params: &UniversalParams<E>,
        poly: &impl MultilinearExtension<Fr>,
        rng: &mut R,
    ) {
        let nv = poly.num_vars();
        assert_ne!(nv, 0);
        let (ck, vk) = MultilinearPC::<E>::trim(&uni_params, nv);
        let point: Vec<_> = (0..nv).map(|_| Fr::rand(rng)).collect();
        let com = MultilinearPC::commit(&ck, poly);
        let proof = MultilinearPC::open(&ck, poly, &point);

        let value = poly.evaluate(&point).unwrap();
        let result = MultilinearPC::check(&vk, &com, &point, value, &proof);
        assert!(result);
    }

    #[test]
    fn setup_commit_verify_correct_polynomials() {
        let mut rng = test_rng();

        // normal polynomials
        let uni_params = MultilinearPC::setup(10, &mut rng);

        let poly1 = DenseMultilinearExtension::rand(8, &mut rng);
        test_polynomial(&uni_params, &poly1, &mut rng);

        let poly2 = SparseMultilinearExtension::rand_with_config(9, 1 << 5, &mut rng);
        test_polynomial(&uni_params, &poly2, &mut rng);

        // single-variate polynomials

        let poly3 = DenseMultilinearExtension::rand(1, &mut rng);
        test_polynomial(&uni_params, &poly3, &mut rng);

        let poly4 = SparseMultilinearExtension::rand_with_config(1, 1 << 1, &mut rng);
        test_polynomial(&uni_params, &poly4, &mut rng);
    }

    #[test]
    #[should_panic]
    fn setup_commit_verify_constant_polynomial() {
        let mut rng = test_rng();

        // normal polynomials
        MultilinearPC::<E>::setup(0, &mut rng);
    }

    #[test]
    fn setup_commit_verify_incorrect_polynomial_should_return_false() {
        let mut rng = test_rng();
        let nv = 8;
        let uni_params = MultilinearPC::setup(nv, &mut rng);
        let poly = DenseMultilinearExtension::rand(nv, &mut rng);
        let nv = uni_params.num_vars;
        let (ck, vk) = MultilinearPC::<E>::trim(&uni_params, nv);
        let point: Vec<_> = (0..nv).map(|_| Fr::rand(&mut rng)).collect();
        let com = MultilinearPC::commit(&ck, &poly);
        let proof = MultilinearPC::open(&ck, &poly, &point);

        let value = poly.evaluate(&point).unwrap();
        let result = MultilinearPC::check(&vk, &com, &point, value + &(1u16.into()), &proof);
        assert!(!result);
    }
}
