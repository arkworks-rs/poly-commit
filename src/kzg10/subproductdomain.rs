use crate::Error;
use ark_ec::PairingEngine;
use ark_ff::{FftField, Field, Zero};
use ark_poly::UVPolynomial;

use ark_poly::polynomial::univariate::DensePolynomial as Poly;

/// Computes the inverse of f mod x^l
pub fn inverse_mod_xl<F: FftField>(f: &Poly<F>, l: usize) -> Option<Poly<F>> {
    let r = ark_std::log2(l);
    let mut g = Poly::<F> {
        coeffs: vec![f.coeffs[0].inverse().unwrap()], //todo unwrap
    };
    let mut i = 2usize;
    for _ in 0..r {
        g = &(&g + &g) - &(f * &(&g * &g)); //todo: g*2?
        g.coeffs.resize(i, F::zero());
        i *= 2;
    }
    Some(g)
}

/// Computes the rev_m(f) function in place
pub fn rev<F: FftField>(f: &mut Poly<F>, m: usize) {
    assert!(f.coeffs.len() - 1 <= m);
    for _ in 0..(m - (f.coeffs.len() - 1)) {
        f.coeffs.push(F::zero());
    }
    f.reverse();
}

/// Divide f by g in nearly linear time
pub fn fast_divide_monic<F: FftField>(f: &Poly<F>, g: &Poly<F>) -> (Poly<F>, Poly<F>) {
    if f.coeffs().len() < g.coeffs().len() {
        return (
            Poly::<F> {
                coeffs: vec![F::zero()],
            },
            f.clone(),
        );
    }
    let m = f.coeffs().len() - g.coeffs().len();

    let mut rev_f = f.clone();
    let mut rev_g = g.clone();
    rev_f.reverse();
    rev_g.reverse();

    let mut q = &rev_f * &inverse_mod_xl::<F>(&rev_g, m + 1).unwrap();
    q.coeffs.resize(m + 1, F::zero());
    rev::<F>(&mut q, m);
    let r = f - &(g * &q);
    (q, r)
}

/// The subproduct tree of a polynomial m over a domain u
#[derive(Debug, Clone)]
pub struct SubproductDomain<F: FftField> {
    /// Domain values
    pub u: Vec<F>,
    /// Subproduct tree over domain u
    pub t: SubproductTree<F>,
    /// Derivative of the subproduct polynomial
    pub prime: Poly<F>, // Derivative
}

impl<F: FftField> SubproductDomain<F> {
    /// Create a subproduct tree domain
    pub fn new(u: Vec<F>) -> SubproductDomain<F> {
        let t = SubproductTree::new(&u);
        let prime = derivative::<F>(&t.m);
        SubproductDomain { u, t, prime }
    }
    /// evaluate a polynomial over the domain
    pub fn evaluate(&self, f: &Poly<F>) -> Vec<F> {
        let mut evals = vec![F::zero(); self.u.len()];
        self.t.evaluate(f, &self.u, &mut evals);
        evals
    }
    /// interpolate a polynomial over the domain
    pub fn interpolate(&self, v: &[F]) -> Poly<F> {
        self.t.interpolate(&self.u, v)
    }
    /// compute the inverse of the lagrange coefficients fast
    pub fn inverse_lagrange_coefficients(&self) -> Vec<F> {
        self.t.inverse_lagrange_coefficients(&self.u)
    }
    /// compute a linear coefficient of lagrange factors times c_i
    pub fn linear_combine(&self, c: &[F]) -> Poly<F> {
        self.t.linear_combine(&self.u, &c)
    }
}

/// A subproduct tree of the subproduct domain
#[derive(Debug, Clone)]
pub struct SubproductTree<F: FftField> {
    /// The left child
    pub left: Option<Box<SubproductTree<F>>>,
    /// The right child
    pub right: Option<Box<SubproductTree<F>>>,
    /// The polynomial for this subdomain
    pub m: Poly<F>,
}

impl<F: FftField> SubproductTree<F> {
    /// Compute the subproduct tree of m = (x - u_0)*...*(x-u_{n-1})
    pub fn new(u: &[F]) -> SubproductTree<F> {
        if u.len() == 1 {
            SubproductTree {
                left: None,
                right: None,
                m: Poly::<F> {
                    coeffs: vec![-u[0], F::one()],
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
    pub fn evaluate(&self, f: &Poly<F>, u: &[F], t: &mut [F]) {
        //todo: assert degree < u.len()
        if u.len() == 1 {
            t[0] = f.coeffs[0];
            return;
        }

        let left = self.left.as_ref().unwrap();
        let right = self.right.as_ref().unwrap();

        let (_, r_0) = fast_divide_monic::<F>(f, &left.m);
        let (_, r_1) = fast_divide_monic::<F>(f, &right.m);

        let n = u.len() / 2;
        let (u_0, u_1) = u.split_at(n);
        let (t_0, t_1) = t.split_at_mut(n);

        left.evaluate(&r_0, u_0, t_0);
        right.evaluate(&r_1, u_1, t_1);
    }
    /// Fast interpolate over this subproduct tree
    pub fn interpolate(&self, u: &[F], v: &[F]) -> Poly<F> {
        let mut lagrange_coeff = self.inverse_lagrange_coefficients(u);

        for (s_i, v_i) in lagrange_coeff.iter_mut().zip(v.iter()) {
            *s_i = s_i.inverse().unwrap() * *v_i;
        }

        self.linear_combine(u, &lagrange_coeff)
    }
    /// Fast compute lagrange coefficients over this subproduct tree
    pub fn inverse_lagrange_coefficients(&self, u: &[F]) -> Vec<F> {
        //assert u.len() == degree of s.m
        if u.len() == 1 {
            return vec![F::one()];
        }
        let mut evals = vec![F::zero(); u.len()];
        let m_prime = derivative::<F>(&self.m);
        self.evaluate(&m_prime, u, &mut evals);
        evals
    }
    /// Fast linear combine over this subproduct tree
    pub fn linear_combine(&self, u: &[F], c: &[F]) -> Poly<F> {
        if u.len() == 1 {
            return Poly::<F> { coeffs: vec![c[0]] };
        }
        let n = u.len() / 2;
        let (u_0, u_1) = u.split_at(n);
        let (c_0, c_1) = c.split_at(n);

        let left = self.left.as_ref().unwrap();
        let right = self.right.as_ref().unwrap();
        let r_0 = left.linear_combine(u_0, c_0);
        let r_1 = right.linear_combine(u_1, c_1);

        &(&right.m * &r_0) + &(&left.m * &r_1)
    }
}
/// compute the derivative of polynomial p
pub fn derivative<F: FftField>(p: &Poly<F>) -> Poly<F> {
    let mut coeffs = Vec::with_capacity(p.coeffs().len() - 1);
    for (i, c) in p.coeffs.iter().enumerate().skip(1) {
        coeffs.push(F::from(i as u64) * c);
    }
    Poly::<F> { coeffs }
}

/// Build a vector representation of the circulant matrix of polynomial p
pub fn build_circulant<F: FftField>(polynomial: &Poly<F>, size: usize) -> Vec<F> {
    let mut circulant = vec![F::zero(); 2 * size];
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
    polynomial: &Poly<E::Fr>,
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
    let mut circulant = build_circulant::<E::Fr>(polynomial, size);

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

    type Fr = <ark_bls12_381::Bls12_381 as PairingEngine>::Fr;
    #[test]
    fn test_inverse() {
        let rng = &mut ark_std::test_rng();

        let degree = 100;
        let l = 101;
        for _ in 0..100 {
            let p = DensePolynomial::<Fr>::rand(degree, rng);
            let p_inv = inverse_mod_xl::<Fr>(&p, l).unwrap();
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
        for g_deg in 1..100 {
            let f = DensePolynomial::<Fr>::rand(degree, rng);
            let mut g = DensePolynomial::<Fr>::rand(g_deg, rng);
            *g.last_mut().unwrap() = Fr::one(); //monic

            let (q, r) = fast_divide_monic::<Fr>(&f, &g);

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

            let s = SubproductDomain::<Fr>::new(points);
            let p = s.interpolate(&evals);

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
            let s = SubproductDomain::<Fr>::new(u);
            let f = s.linear_combine(&c);

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
            let s = SubproductDomain::<Fr>::new(u);
            let f = s.inverse_lagrange_coefficients();

            for (i, j) in s.u.iter().zip(f.iter()) {
                assert_eq!(s.prime.evaluate(i), *j);
            }
        }
    }
}
