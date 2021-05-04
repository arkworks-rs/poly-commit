use crate::Error;
use ark_ec::PairingEngine;
use ark_ff::{Field, One, Zero};
use ark_poly::UVPolynomial;

use ark_poly::polynomial::univariate::DensePolynomial;

type P<E: PairingEngine> = DensePolynomial<E::Fr>;

/// Computes the inverse of f mod x^l
pub fn inverse_mod_xl<E: PairingEngine>(f: &P<E>, l: usize) -> Option<P<E>> {
    //use std::ops::Mul;
    //let r =
    //    std::mem::size_of::<u64>() * 8 - (l as u64).leading_zeros() as usize; // ceil(log_2(l))

    //assert_eq!((l as f64).log2().ceil() as usize, r);
    let r = (l as f64).log2().ceil() as usize; //TODO: rounding problems??
    let mut g = DensePolynomial::<E::Fr> {
        coeffs: vec![f.coeffs[0].inverse().unwrap()], //todo unwrap
    };
    let mut i = 2usize;
    for _ in 0..r {
        g = &(&g + &g) - &(f * &(&g * &g)); //todo: g*2?
        g.coeffs.resize(i, E::Fr::zero());
        i *= 2;
    }
    Some(g)
}

/// Computes the rev_m(f) function in place
pub fn rev<E: PairingEngine>(f: &mut P<E>, m: usize) {
    assert!(f.coeffs.len() - 1 <= m);
    for _ in 0..(m - (f.coeffs.len() - 1)) {
        f.coeffs.push(E::Fr::zero());
    }
    f.reverse();
}

/// Divide f by g in nearly linear time
pub fn fast_divide_monic<E: PairingEngine>(f: &P<E>, g: &P<E>) -> (P<E>, P<E>) {
    if f.coeffs().len() < g.coeffs().len() {
        return (
            P::<E> {
                coeffs: vec![E::Fr::zero()],
            },
            f.clone(),
        );
    }
    let m = f.coeffs().len() - g.coeffs().len();

    let mut rev_f = f.clone();
    let mut rev_g = g.clone();
    rev_f.reverse();
    rev_g.reverse();

    let mut q = &rev_f * &inverse_mod_xl::<E>(&rev_g, m + 1).unwrap();
    q.coeffs.resize(m + 1, E::Fr::zero());
    rev::<E>(&mut q, m);
    let r = f - &(g * &q);
    (q, r)
}

/// The subproduct tree of a polynomial m over a domain u
#[derive(Debug, Clone)]
pub struct SubproductDomain<E: PairingEngine> {
    /// Domain values
    pub u: Vec<E::Fr>,
    /// Subproduct tree over domain u
    pub t: SubproductTree<E>,
    /// Derivative of the subproduct polynomial
    pub prime: P<E>, // Derivative
}

impl<E: PairingEngine> SubproductDomain<E> {
    /// Create a subproduct tree domain
    pub fn new(u: Vec<E::Fr>) -> SubproductDomain<E> {
        let t = SubproductTree::new(&u);
        let prime = derivative::<E>(&t.m);
        SubproductDomain { u, t, prime }
    }
    /// evaluate a polynomial over the domain
    pub fn fast_evaluate(&self, f: &P<E>) -> Vec<E::Fr> {
        let mut evals = vec![E::Fr::zero(); self.u.len()];
        self.t.fast_evaluate(f, &self.u, &mut evals);
        evals
    }
    /// interpolate a polynomial over the domain
    pub fn fast_interpolate(&self, v: &[E::Fr]) -> P<E> {
        self.t.fast_interpolate(&self.u, v)
    }
    /// compute the inverse of the lagrange coefficients fast
    pub fn fast_inverse_lagrange_coefficients(&self) -> Vec<E::Fr> {
        self.t.fast_inverse_lagrange_coefficients(&self.u)
    }
    /// compute a linear coefficient of lagrange factors times c_i
    pub fn fast_linear_combine(&self, c: &[E::Fr]) -> P<E> {
        self.t.fast_linear_combine(&self.u, &c)
    }
}

/// A subproduct tree of the subproduct domain
#[derive(Debug, Clone)]
pub struct SubproductTree<E: PairingEngine> {
    /// The left child
    pub left: Option<Box<SubproductTree<E>>>,
    /// The right child
    pub right: Option<Box<SubproductTree<E>>>,
    /// The polynomial for this subdomain
    pub m: P<E>,
}

impl<E: PairingEngine> SubproductTree<E> {
    /// Compute the subproduct tree of m = (x - u_0)*...*(x-u_{n-1})
    pub fn new(u: &[E::Fr]) -> SubproductTree<E> {
        if u.len() == 1 {
            SubproductTree {
                left: None,
                right: None,
                m: P::<E> {
                    coeffs: vec![-u[0], E::Fr::one()],
                },
            }
        } else {
            let n = u.len() / 2;
            let (u_0, u_1) = u.split_at(n);
            let left = Box::new(SubproductTree::new(u_0));
            let right = Box::new(SubproductTree::new(u_1));
            let m = &left.m * &right.m;
            SubproductTree {
                left: Some(left),
                right: Some(right),
                m,
            }
        }
    }
    /// Fast evaluate f over this subproduct tree
    pub fn fast_evaluate(&self, f: &P<E>, u: &[E::Fr], t: &mut [E::Fr]) {
        //todo: assert degree < u.len()
        if u.len() == 1 {
            t[0] = f.coeffs[0];
            return;
        }

        let left = self.left.as_ref().unwrap();
        let right = self.right.as_ref().unwrap();

        let (q_0, r_0) = fast_divide_monic::<E>(f, &left.m);
        let (_, r_1) = fast_divide_monic::<E>(f, &right.m);

        let n = u.len() / 2;
        let (u_0, u_1) = u.split_at(n);
        let (t_0, t_1) = t.split_at_mut(n);

        left.fast_evaluate(&r_0, u_0, t_0);
        right.fast_evaluate(&r_1, u_1, t_1);
    }
    /// Fast interpolate over this subproduct tree
    pub fn fast_interpolate(&self, u: &[E::Fr], v: &[E::Fr]) -> P<E> {
        let mut lagrange_coeff = self.fast_inverse_lagrange_coefficients(u);

        for (s_i, v_i) in lagrange_coeff.iter_mut().zip(v.iter()) {
            *s_i = s_i.inverse().unwrap() * *v_i;
        }

        self.fast_linear_combine(u, &lagrange_coeff)
    }
    /// Fast compute lagrange coefficients over this subproduct tree
    pub fn fast_inverse_lagrange_coefficients(&self, u: &[E::Fr]) -> Vec<E::Fr> {
        //assert u.len() == degree of s.m
        if u.len() == 1 {
            return vec![E::Fr::one()];
        }
        let mut evals = vec![E::Fr::zero(); u.len()];
        let m_prime = derivative::<E>(&self.m);
        self.fast_evaluate(&m_prime, u, &mut evals);
        evals
    }
    /// Fast linear combine over this subproduct tree
    pub fn fast_linear_combine(&self, u: &[E::Fr], c: &[E::Fr]) -> P<E> {
        if u.len() == 1 {
            return P::<E> { coeffs: vec![c[0]] };
        }
        let n = u.len() / 2;
        let (u_0, u_1) = u.split_at(n);
        let (c_0, c_1) = c.split_at(n);

        let left = self.left.as_ref().unwrap();
        let right = self.right.as_ref().unwrap();
        let r_0 = left.fast_linear_combine(u_0, c_0);
        let r_1 = right.fast_linear_combine(u_1, c_1);

        &(&right.m * &r_0) + &(&left.m * &r_1)
    }
}
/// compute the derivative of polynomial p
pub fn derivative<E: PairingEngine>(p: &P<E>) -> P<E> {
    let mut coeffs = Vec::with_capacity(p.coeffs().len() - 1);
    for (i, c) in p.coeffs.iter().enumerate().skip(1) {
        coeffs.push(E::Fr::from(i as u64) * c);
    }
    P::<E> { coeffs }
}

/// Build a vector representation of the circulant matrix of polynomial p
pub fn build_circulant<E: PairingEngine>(polynomial: &P<E>, size: usize) -> Vec<E::Fr> {
    let mut circulant = vec![E::Fr::zero(); 2 * size];
    let coeffs = polynomial.coeffs();
    if size == coeffs.len() - 1 {
        circulant[0] = *coeffs.last().unwrap();
        circulant[size] = *coeffs.last().unwrap();
        circulant[size + 1..size + 1 + coeffs.len() - 2]
            .copy_from_slice(&coeffs[1..coeffs.len() - 1]);
    } else {
        circulant[size + 1..size + 1 + coeffs.len() - 1].copy_from_slice(&coeffs[1..]);
    }
    circulant
}
/// Computes the Toeplitz matrix of polynomial times the vector v
pub fn toeplitz_mul<E: PairingEngine>(
    polynomial: &P<E>,
    v: &[E::G1Affine],
    size: usize,
) -> Result<(Vec<E::G1Projective>, E::Fr), Error> {
    use ark_ec::AffineCurve;
    use ark_poly::EvaluationDomain;
    let m = polynomial.coeffs().len() - 1;
    let size = ark_std::cmp::max(size, m);

    let domain = ark_poly::Radix2EvaluationDomain::<E::Fr>::new(2 * size)
        .ok_or(Error::AmortizedOpeningTooLarge(size))?;

    let size = domain.size() / 2;
    let mut circulant = build_circulant::<E>(polynomial, size);

    let mut tmp: Vec<E::G1Projective> = Vec::with_capacity(domain.size());

    for _ in 0..(size - v.len()) {
        tmp.push(E::G1Projective::zero());
    }

    for i in v.iter().rev() {
        tmp.push(i.into_projective());
    }

    tmp.resize(domain.size(), E::G1Projective::zero());
    domain.fft_in_place(&mut tmp);
    domain.fft_in_place(&mut circulant);

    for (i, j) in tmp.iter_mut().zip(circulant.iter()) {
        *i *= *j;
    }

    // NOTE: it is possible to avoid scaling by domain_size_inv, and defer this to later
    domain.ifft_in_place(&mut tmp);

    Ok((
        tmp[..size].to_vec(),
        E::Fr::from(domain.size() as u64).inverse().unwrap(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ec::PairingEngine;
    use ark_ff::{One, Zero};
    use ark_poly::polynomial::univariate::DensePolynomial;
    use ark_poly::Polynomial;
    use ark_poly::UVPolynomial;
    use ark_std::UniformRand;

    type E = ark_bls12_381::Bls12_381;
    type Fr = <ark_bls12_381::Bls12_381 as PairingEngine>::Fr;
    #[test]
    fn test_inverse() {
        let rng = &mut ark_std::test_rng();

        let degree = 100;
        let l = 101;
        for _ in 0..100 {
            let p = DensePolynomial::<Fr>::rand(degree, rng);
            let p_inv = inverse_mod_xl::<E>(&p, l).unwrap();
            let mut t = &p * &p_inv;
            t.coeffs.resize(l, Fr::zero());
            assert_eq!(t.coeffs[0], Fr::one());
            for i in t.iter().skip(1) {
                assert_eq!(*i, Fr::zero());
            }
        }
    }

    #[test]
    fn test_divide() {
        let rng = &mut ark_std::test_rng();

        let degree = 100;
        let l = 101;
        for g_deg in 1..100 {
            let f = DensePolynomial::<Fr>::rand(degree, rng);
            let mut g = DensePolynomial::<Fr>::rand(g_deg, rng);
            *g.last_mut().unwrap() = Fr::one(); //monic

            let (q, r) = fast_divide_monic::<E>(&f, &g);

            let t = &(&q * &g) + &r;

            for (i, j) in t.coeffs.iter().zip(f.coeffs.iter()) {
                assert_eq!(*i, *j);
            }
        }
    }

    #[test]
    fn test_interpolate() {
        let rng = &mut ark_std::test_rng();
        for d in 1..100 {
            let mut points = vec![];
            let mut evals = vec![];
            for _ in 0..d {
                points.push(Fr::rand(rng));
                evals.push(Fr::rand(rng));
            }

            let s = SubproductDomain::<E>::new(points);
            let p = s.fast_interpolate(&evals);

            for (x, y) in s.u.iter().zip(evals.iter()) {
                assert_eq!(p.evaluate(x), *y)
            }
        }
    }

    #[test]
    fn test_linear_combine() {
        let rng = &mut ark_std::test_rng();
        for d in 1..100 {
            let mut u = vec![];
            let mut c = vec![];
            for _ in 0..d {
                u.push(Fr::rand(rng));
                c.push(Fr::rand(rng));
            }
            let s = SubproductDomain::<E>::new(u);
            let f = s.fast_linear_combine(&c);

            let r = Fr::rand(rng);
            let m = s.t.m.evaluate(&r);
            let mut total = Fr::zero();
            for (u_i, c_i) in s.u.iter().zip(c.iter()) {
                total += m * *c_i / (r - u_i);
            }
            assert_eq!(f.evaluate(&r), total);
        }
    }

    #[test]
    fn test_inv_lagrange() {
        let rng = &mut ark_std::test_rng();
        for d in 1..100 {
            let mut u = vec![];
            for _ in 0..d {
                u.push(Fr::rand(rng));
            }
            let s = SubproductDomain::<E>::new(u);
            let f = s.fast_inverse_lagrange_coefficients();

            for (a, (i, j)) in s.u.iter().zip(f.iter()).enumerate() {
                assert_eq!(s.prime.evaluate(i), *j);
            }
        }
    }
}
