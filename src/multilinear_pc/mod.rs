use ark_ec::msm::{FixedBaseMSM, VariableBaseMSM};
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{Field, PrimeField};
use ark_ff::{One, Zero};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::collections::LinkedList;
use ark_std::iter::FromIterator;
use ark_std::marker::PhantomData;
use ark_std::UniformRand;
use rand_core::RngCore;

use crate::multilinear_pc::data_structures::{
    Commitment, CommitterKey, Proof, UniversalParams, VerifierKey,
};
use crate::Error;

/// data structures used by multilinear extension commitment scheme
pub mod data_structures;

/// Polynomial Commitment Scheme on multilinear extensions.
pub struct MultilinearPC<E: PairingEngine, P: MultilinearExtension<E::Fr>> {
    _engine: PhantomData<E>,
    _polynomial: PhantomData<P>,
}

impl<E: PairingEngine, P: MultilinearExtension<E::Fr>> MultilinearPC<E, P> {
    /// setup
    pub fn setup<R: RngCore>(num_vars: usize, rng: &mut R) -> Result<UniversalParams<E>, Error> {
        let g: E::G1Projective = E::G1Projective::rand(rng);
        let h: E::G2Projective = E::G2Projective::rand(rng);
        let g = g.into_affine();
        let h = h.into_affine();
        let mut powers_of_g = Vec::new();
        let mut powers_of_h = Vec::new();
        let t: Vec<_> = (0..num_vars).map(|_| E::Fr::rand(rng)).collect();
        let scalar_bits = E::Fr::size_in_bits();

        let mut eq: LinkedList<DenseMultilinearExtension<E::Fr>> =
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
        let window_size = FixedBaseMSM::get_mul_window_size(total_scalars);
        let g_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, g.into_projective());
        let h_table = FixedBaseMSM::get_window_table(scalar_bits, window_size, h.into_projective());

        let pp_g = E::G1Projective::batch_normalization_into_affine(
            &FixedBaseMSM::multi_scalar_mul(scalar_bits, window_size, &g_table, &pp_powers),
        );
        let pp_h = E::G2Projective::batch_normalization_into_affine(
            &FixedBaseMSM::multi_scalar_mul(scalar_bits, window_size, &h_table, &pp_powers),
        );
        let mut start = 0;
        for i in 0..num_vars {
            let size = 1 << (num_vars - i);
            let pp_k_g = (&pp_g[start..(start + size)]).to_vec();
            let pp_k_h = (&pp_h[start..(start + size)]).to_vec();
            powers_of_g.push(pp_k_g);
            powers_of_h.push(pp_k_h);
            start += size;
        }

        // end_timer!(variable_mul_timer);
        // calculate vp
        // let vp_generation_timer = start_timer!(|| "VP generation");
        let g_mask = {
            let window_size = FixedBaseMSM::get_mul_window_size(num_vars);
            let g_table =
                FixedBaseMSM::get_window_table(scalar_bits, window_size, g.into_projective());
            E::G1Projective::batch_normalization_into_affine(&FixedBaseMSM::multi_scalar_mul(
                scalar_bits,
                window_size,
                &g_table,
                &t,
            ))
        };
        // end_timer!(vp_generation_timer);

        Ok(UniversalParams {
            nv: num_vars,
            g,
            g_mask,
            h,
            powers_of_g,
            powers_of_h,
        })
    }

    /// get committer key and verifier key from universal parameters
    pub fn get_keys(params: &UniversalParams<E>) -> (CommitterKey<E>, VerifierKey<E>) {
        let ck = CommitterKey {
            powers_of_h: params.powers_of_h.to_vec(),
            powers_of_g: params.powers_of_g.to_vec(),
            g: params.g.clone(),
            h: params.h.clone(),
            nv: params.nv,
        };
        let vk = VerifierKey {
            nv: params.nv,
            g: params.g.clone(),
            h: params.h.clone(),
            g_mask_random: params.g_mask.clone(),
        };
        (ck, vk)
    }

    /// commit
    pub fn commit(
        ck: &CommitterKey<E>,
        polynomial: &impl MultilinearExtension<E::Fr>,
    ) -> Commitment<E> {
        let nv = polynomial.num_vars();
        let scalars: Vec<_> = polynomial
            .to_evaluations()
            .into_iter()
            .map(|x| x.into_repr())
            .collect();
        let g_product =
            VariableBaseMSM::multi_scalar_mul(&ck.powers_of_g[0], scalars.as_slice()).into_affine();
        Commitment { nv, g_product }
    }

    /// On input a polynomial `p` and a point `point`, outputs a proof for the same.
    pub fn open(
        ck: &CommitterKey<E>,
        polynomial: &impl MultilinearExtension<E::Fr>,
        point: &[E::Fr],
    ) -> Proof<E> {
        assert_eq!(polynomial.num_vars(), ck.nv, "Invalid size of polynomial");
        let nv = polynomial.num_vars();
        let mut r: Vec<Vec<E::Fr>> = (0..nv + 1).map(|_| Vec::new()).collect();
        let mut q: Vec<Vec<E::Fr>> = (0..nv + 1).map(|_| Vec::new()).collect();

        r[nv] = polynomial.to_evaluations();

        let mut proofs = Vec::new();
        for i in 0..nv {
            let k = nv - i;
            let point_at_k = point[i];
            q[k] = (0..(1 << (k - 1))).map(|_| E::Fr::zero()).collect();
            r[k - 1] = (0..(1 << (k - 1))).map(|_| E::Fr::zero()).collect();
            for b in 0..(1 << (k - 1)) {
                q[k][b] = r[k][(b << 1) + 1] - &r[k][b << 1];
                r[k - 1][b] = r[k][b << 1] * &(E::Fr::one() - &point_at_k)
                    + &(r[k][(b << 1) + 1] * &point_at_k);
            }
            let scalars: Vec<_> = (0..(1 << k))
                .map(|x| q[k][x >> 1].into_repr()) // fine
                .collect();

            let pi_h =
                VariableBaseMSM::multi_scalar_mul(&ck.powers_of_h[i], &scalars).into_affine(); // no need to move outside and partition
            proofs.push(pi_h);
        }

        Proof { proofs }
    }

    /// Verifies that `value` is the evaluation at `x` of the polynomial
    /// committed inside `comm`.
    pub fn check<'a>(
        vk: &VerifierKey<E>,
        commitment: &Commitment<E>,
        point: &[E::Fr],
        value: E::Fr,
        proof: &Proof<E>,
    ) -> bool {
        let left = E::pairing(
            commitment.g_product.into_projective() - &vk.g.mul(value),
            vk.h,
        );

        let scalar_size = E::Fr::size_in_bits();
        let window_size = FixedBaseMSM::get_mul_window_size(vk.nv);
        // let timer = start_timer!(|| "MSM");
        let g_table =
            FixedBaseMSM::get_window_table(scalar_size, window_size, vk.g.into_projective());
        let g_mul: Vec<E::G1Projective> =
            FixedBaseMSM::multi_scalar_mul(scalar_size, window_size, &g_table, point);

        let pairing_lefts: Vec<_> = (0..vk.nv)
            .map(|i| vk.g_mask_random[i].into_projective() - &g_mul[i])
            .collect();
        let pairing_lefts: Vec<E::G1Affine> =
            E::G1Projective::batch_normalization_into_affine(&pairing_lefts);
        let pairing_lefts: Vec<E::G1Prepared> = pairing_lefts
            .into_iter()
            .map(|x| E::G1Prepared::from(x))
            .collect();

        let pairing_rights: Vec<E::G2Prepared> = proof
            .proofs
            .iter()
            .map(|x| E::G2Prepared::from(*x))
            .collect();

        let pairings: Vec<_> = pairing_lefts
            .into_iter()
            .zip(pairing_rights.into_iter())
            .collect();
        let right = E::product_of_pairings(pairings.iter());
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
