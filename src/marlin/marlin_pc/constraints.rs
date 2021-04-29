use crate::{
    constraints::{EvaluationsVar, LabeledPointVar, PCCheckRandomDataVar, PCCheckVar, QuerySetVar},
    data_structures::LabeledCommitment,
    kzg10::{Proof, VerifierKey as KZG10VerifierKey},
    marlin_pc::{
        data_structures::{Commitment, VerifierKey},
        MarlinKZG10, PreparedCommitment, PreparedVerifierKey,
    },
    BTreeMap, BTreeSet, BatchLCProof, LinearCombinationCoeffVar, LinearCombinationVar,
    PrepareGadget, String, ToString, Vec,
};
use ark_ec::{AffineCurve, CurveCycle, PairingEngine, PairingFriendlyCycle};
use ark_ff::{PrimeField, ToConstraintField};
use ark_nonnative_field::{NonNativeFieldMulResultVar, NonNativeFieldVar};
use ark_poly::UVPolynomial;
use ark_r1cs_std::{fields::fp::FpVar, prelude::*, ToConstraintFieldGadget};
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, Result as R1CSResult, SynthesisError};
use ark_std::{borrow::Borrow, convert::TryInto, marker::PhantomData, ops::Div, vec};

type E2Fr<E> = <<E as CurveCycle>::E2 as AffineCurve>::ScalarField;
type E2Fq<E> = <<E as CurveCycle>::E2 as AffineCurve>::BaseField;

/// High level variable representing the verification key of the `MarlinKZG10` polynomial commitment
/// scheme.
#[derive(Derivative)]
#[derivative(Clone(bound = "E: PairingFriendlyCycle"))]
pub struct VerifierKeyVar<E: PairingFriendlyCycle, PG: PairingVar<E::Engine2, E2Fq<E>>>
where
    E2Fq<E>: PrimeField,
{
    /// Generator of G1.
    pub g: PG::G1Var,
    /// Generator of G2.
    pub h: PG::G2Var,
    /// Generator of G1, times first monomial.
    pub beta_h: PG::G2Var,
    /// Used for the shift powers associated with different degree bounds.
    pub degree_bounds_and_shift_powers: Option<Vec<(usize, FpVar<E2Fq<E>>, PG::G1Var)>>,
}

impl<E, PG> VerifierKeyVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    /// Find the appropriate shift for the degree bound.
    #[tracing::instrument(target = "r1cs", skip(self, cs))]
    pub fn get_shift_power(
        &self,
        cs: impl Into<Namespace<E2Fq<E>>>,
        bound: &FpVar<E2Fq<E>>,
    ) -> Option<PG::G1Var> {
        // Search the bound using PIR
        if self.degree_bounds_and_shift_powers.is_none() {
            None
        } else {
            let ns = cs.into();
            let cs = ns.cs();

            let degree_bounds_and_shift_powers =
                self.degree_bounds_and_shift_powers.clone().unwrap();

            let mut pir_vector = vec![false; degree_bounds_and_shift_powers.len()];

            let desired_bound_value = bound.value().unwrap_or_default();

            for (i, (_, this_bound, _)) in degree_bounds_and_shift_powers.iter().enumerate() {
                if this_bound
                    .value()
                    .unwrap_or_default()
                    .eq(&desired_bound_value)
                {
                    pir_vector[i] = true;
                    break;
                }
            }

            let mut pir_vector_gadgets = Vec::new();
            for bit in pir_vector.iter() {
                pir_vector_gadgets.push(
                    Boolean::new_witness(ark_relations::ns!(cs, "alloc_pir"), || Ok(bit)).unwrap(),
                );
            }

            // Sum of the PIR values are equal to one
            let mut sum = FpVar::<E2Fq<E>>::zero();
            let one = FpVar::<E2Fq<E>>::one();
            for pir_gadget in pir_vector_gadgets.iter() {
                sum += &FpVar::<E2Fq<E>>::from(pir_gadget.clone());
            }
            sum.enforce_equal(&one).unwrap();

            // PIR the value
            let mut found_bound = FpVar::<E2Fq<E>>::zero();

            let mut found_shift_power = PG::G1Var::zero();

            for (pir_gadget, (_, degree, shift_power)) in pir_vector_gadgets
                .iter()
                .zip(degree_bounds_and_shift_powers.iter())
            {
                found_bound =
                    FpVar::<E2Fq<E>>::conditionally_select(pir_gadget, degree, &found_bound)
                        .unwrap();

                found_shift_power =
                    PG::G1Var::conditionally_select(pir_gadget, shift_power, &found_shift_power)
                        .unwrap();
            }

            found_bound.enforce_equal(&bound).unwrap();

            Some(found_shift_power)
        }
    }
}

impl<E, PG> AllocVar<VerifierKey<E::Engine2>, E2Fq<E>> for VerifierKeyVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, val))]
    fn new_variable<T>(
        cs: impl Into<Namespace<E2Fq<E>>>,
        val: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> R1CSResult<Self>
    where
        T: Borrow<VerifierKey<E::Engine2>>,
    {
        let vk_orig = val()?.borrow().clone();

        let ns = cs.into();
        let cs = ns.cs();

        let VerifierKey {
            vk,
            degree_bounds_and_shift_powers,
            ..
        } = vk_orig;

        let degree_bounds_and_shift_powers = degree_bounds_and_shift_powers.map(|vec| {
            vec.iter()
                .map(|(s, g)| {
                    (
                        *s,
                        FpVar::<E2Fq<E>>::new_variable(
                            ark_relations::ns!(cs, "degree bound"),
                            || Ok(<E2Fq<E> as From<u128>>::from(*s as u128)),
                            mode,
                        )
                        .unwrap(),
                        PG::G1Var::new_variable(ark_relations::ns!(cs, "pow"), || Ok(*g), mode)
                            .unwrap(),
                    )
                })
                .collect()
        });

        let KZG10VerifierKey { g, h, beta_h, .. } = vk;

        let g = PG::G1Var::new_variable(ark_relations::ns!(cs, "g"), || Ok(g), mode)?;
        let h = PG::G2Var::new_variable(ark_relations::ns!(cs, "h"), || Ok(h), mode)?;
        let beta_h =
            PG::G2Var::new_variable(ark_relations::ns!(cs, "beta_h"), || Ok(beta_h), mode)?;

        Ok(Self {
            g,
            h,
            beta_h,
            degree_bounds_and_shift_powers,
        })
    }
}

impl<E, PG> ToBytesGadget<E2Fq<E>> for VerifierKeyVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(self))]
    fn to_bytes(&self) -> R1CSResult<Vec<UInt8<E2Fq<E>>>> {
        let mut bytes = Vec::new();

        bytes.extend_from_slice(&self.g.to_bytes()?);
        bytes.extend_from_slice(&self.h.to_bytes()?);
        bytes.extend_from_slice(&self.beta_h.to_bytes()?);

        if self.degree_bounds_and_shift_powers.is_some() {
            let degree_bounds_and_shift_powers =
                self.degree_bounds_and_shift_powers.as_ref().unwrap();
            for (_, degree_bound, shift_power) in degree_bounds_and_shift_powers.iter() {
                bytes.extend_from_slice(&degree_bound.to_bytes()?);
                bytes.extend_from_slice(&shift_power.to_bytes()?);
            }
        }

        Ok(bytes)
    }
}

impl<E, PG> ToConstraintFieldGadget<E2Fq<E>> for VerifierKeyVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
    PG::G1Var: ToConstraintFieldGadget<E2Fq<E>>,
    PG::G2Var: ToConstraintFieldGadget<E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(self))]
    fn to_constraint_field(&self) -> R1CSResult<Vec<FpVar<E2Fq<E>>>> {
        let mut res = Vec::new();

        let mut g_gadget = self.g.to_constraint_field()?;
        let mut h_gadget = self.h.to_constraint_field()?;
        let mut beta_h_gadget = self.beta_h.to_constraint_field()?;

        res.append(&mut g_gadget);
        res.append(&mut h_gadget);
        res.append(&mut beta_h_gadget);

        if self.degree_bounds_and_shift_powers.as_ref().is_some() {
            let list = self.degree_bounds_and_shift_powers.as_ref().unwrap();
            for (_, d_gadget, shift_power) in list.iter() {
                let mut d_elems = d_gadget.to_constraint_field()?;
                let mut shift_power_elems = shift_power.to_constraint_field()?;

                res.append(&mut d_elems);
                res.append(&mut shift_power_elems);
            }
        }

        Ok(res)
    }
}

/// High level variable representing the verification key of the `MarlinKZG10` polynomial commitment
/// scheme, prepared for use in arithmetic.
#[allow(clippy::type_complexity)]
#[derive(Derivative)]
#[derivative(Clone(bound = "E: PairingFriendlyCycle"))]
pub struct PreparedVerifierKeyVar<E: PairingFriendlyCycle, PG: PairingVar<E::Engine2, E2Fq<E>>>
where
    E2Fq<E>: PrimeField,
{
    /// Generator of G1.
    pub prepared_g: Vec<PG::G1Var>,
    /// Generator of G2.
    pub prepared_h: PG::G2PreparedVar,
    /// Generator of G1, times first monomial.
    pub prepared_beta_h: PG::G2PreparedVar,
    /// Used for the shift powers associated with different degree bounds.
    pub prepared_degree_bounds_and_shift_powers:
        Option<Vec<(usize, FpVar<E2Fq<E>>, Vec<PG::G1Var>)>>,
    /// Indicate whether or not it is a constant allocation (which decides whether or not shift
    /// powers are precomputed).
    pub constant_allocation: bool,
    /// If not a constant allocation, the original vk is attached (for computing the shift power
    /// series).
    pub origin_vk: Option<VerifierKeyVar<E, PG>>,
}

impl<E, PG> PreparedVerifierKeyVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    /// Find the appropriate shift for the degree bound.
    pub fn get_shift_power(
        &self,
        cs: impl Into<Namespace<E2Fq<E>>>,
        bound: &FpVar<E2Fq<E>>,
    ) -> Option<Vec<PG::G1Var>> {
        if self.constant_allocation {
            if self.prepared_degree_bounds_and_shift_powers.is_none() {
                None
            } else {
                let prepared_degree_bounds_and_shift_powers = self
                    .prepared_degree_bounds_and_shift_powers
                    .as_ref()
                    .unwrap();
                let bound_value = bound.value().unwrap_or_default();

                for (_, bound, prepared_shift_powers) in
                    prepared_degree_bounds_and_shift_powers.iter()
                {
                    if bound.value().unwrap_or_default() == bound_value {
                        return Some((*prepared_shift_powers).clone());
                    }
                }

                None
            }
        } else {
            let shift_power = self.origin_vk.as_ref().unwrap().get_shift_power(cs, bound);

            if let Some(shift_power) = shift_power {
                let mut prepared_shift_gadgets = Vec::<PG::G1Var>::new();

                let supported_bits = E2Fr::<E>::size_in_bits();

                let mut cur: PG::G1Var = shift_power;
                for _ in 0..supported_bits {
                    prepared_shift_gadgets.push(cur.clone());
                    cur.double_in_place().unwrap();
                }

                Some(prepared_shift_gadgets)
            } else {
                None
            }
        }
    }
}

impl<E, PG> PrepareGadget<VerifierKeyVar<E, PG>, E2Fq<E>> for PreparedVerifierKeyVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(unprepared))]
    fn prepare(unprepared: &VerifierKeyVar<E, PG>) -> R1CSResult<Self> {
        let supported_bits = <E2Fr<E> as PrimeField>::size_in_bits();
        let mut prepared_g = Vec::<PG::G1Var>::new();

        let mut g: PG::G1Var = unprepared.g.clone();
        for _ in 0..supported_bits {
            prepared_g.push(g.clone());
            g.double_in_place()?;
        }

        let prepared_h = PG::prepare_g2(&unprepared.h)?;
        let prepared_beta_h = PG::prepare_g2(&unprepared.beta_h)?;

        let prepared_degree_bounds_and_shift_powers =
            if unprepared.degree_bounds_and_shift_powers.is_some() {
                let mut res = Vec::<(usize, FpVar<E2Fq<E>>, Vec<PG::G1Var>)>::new();

                for (d, d_gadget, shift_power) in unprepared
                    .degree_bounds_and_shift_powers
                    .as_ref()
                    .unwrap()
                    .iter()
                {
                    res.push((*d, (*d_gadget).clone(), vec![shift_power.clone()]));
                }

                Some(res)
            } else {
                None
            };

        Ok(Self {
            prepared_g,
            prepared_h,
            prepared_beta_h,
            prepared_degree_bounds_and_shift_powers,
            constant_allocation: false,
            origin_vk: Some(unprepared.clone()),
        })
    }
}

impl<E, PG> AllocVar<PreparedVerifierKey<E::Engine2>, E2Fq<E>> for PreparedVerifierKeyVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable<T>(
        cs: impl Into<Namespace<E2Fq<E>>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> R1CSResult<Self>
    where
        T: Borrow<PreparedVerifierKey<E::Engine2>>,
    {
        let t = f()?;
        let obj = t.borrow();

        let ns = cs.into();
        let cs = ns.cs();

        let mut prepared_g = Vec::<PG::G1Var>::new();
        for g in obj.prepared_vk.prepared_g.iter() {
            prepared_g.push(<PG::G1Var as AllocVar<
                <E::Engine2 as PairingEngine>::G1Affine,
                E2Fq<E>,
            >>::new_variable(
                ark_relations::ns!(cs, "g"), || Ok(*g), mode
            )?);
        }

        let prepared_h = PG::G2PreparedVar::new_variable(
            ark_relations::ns!(cs, "h"),
            || Ok(&obj.prepared_vk.prepared_h),
            mode,
        )?;
        let prepared_beta_h = PG::G2PreparedVar::new_variable(
            ark_relations::ns!(cs, "beta_h"),
            || Ok(&obj.prepared_vk.prepared_beta_h),
            mode,
        )?;

        let prepared_degree_bounds_and_shift_powers =
            if obj.prepared_degree_bounds_and_shift_powers.is_some() {
                let mut res = Vec::<(usize, FpVar<E2Fq<E>>, Vec<PG::G1Var>)>::new();

                for (d, shift_power_elems) in obj
                    .prepared_degree_bounds_and_shift_powers
                    .as_ref()
                    .unwrap()
                    .iter()
                {
                    let mut gadgets = Vec::<PG::G1Var>::new();
                    for shift_power_elem in shift_power_elems.iter() {
                        gadgets.push(<PG::G1Var as AllocVar<
                            <E::Engine2 as PairingEngine>::G1Affine,
                            E2Fq<E>,
                        >>::new_variable(
                            cs.clone(), || Ok(shift_power_elem), mode
                        )?);
                    }

                    let d_gadget = FpVar::<E2Fq<E>>::new_variable(
                        cs.clone(),
                        || Ok(<E2Fq<E> as From<u128>>::from(*d as u128)),
                        mode,
                    )?;

                    res.push((*d, d_gadget, gadgets));
                }
                Some(res)
            } else {
                None
            };

        Ok(Self {
            prepared_g,
            prepared_h,
            prepared_beta_h,
            prepared_degree_bounds_and_shift_powers,
            constant_allocation: true,
            origin_vk: None,
        })
    }
}

/// High level variable representing a commitment in the `MarlinKZG10` polynomial commitment scheme.
#[derive(Derivative)]
#[derivative(Clone(bound = "E: PairingFriendlyCycle"))]
pub struct CommitmentVar<E: PairingFriendlyCycle, PG: PairingVar<E::Engine2, E2Fq<E>>>
where
    E2Fq<E>: PrimeField,
{
    comm: PG::G1Var,
    shifted_comm: Option<PG::G1Var>,
}

impl<E, PG> AllocVar<Commitment<E::Engine2>, E2Fq<E>> for CommitmentVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, value_gen))]
    fn new_variable<T>(
        cs: impl Into<Namespace<E2Fq<E>>>,
        value_gen: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> R1CSResult<Self>
    where
        T: Borrow<Commitment<E::Engine2>>,
    {
        value_gen().and_then(|commitment| {
            let ns = cs.into();
            let cs = ns.cs();

            let commitment = *commitment.borrow();
            let comm = commitment.comm;
            let comm_gadget = PG::G1Var::new_variable(cs.clone(), || Ok(comm.0), mode)?;

            let shifted_comm = commitment.shifted_comm;
            let shifted_comm_gadget = if let Some(shifted_comm) = shifted_comm {
                Some(PG::G1Var::new_variable(cs, || Ok(shifted_comm.0), mode)?)
            } else {
                None
            };

            Ok(Self {
                comm: comm_gadget,
                shifted_comm: shifted_comm_gadget,
            })
        })
    }
}

impl<E, PG> ToConstraintFieldGadget<E2Fq<E>> for CommitmentVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    E::E2: ToConstraintField<E2Fq<E>>,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
    PG::G1Var: ToConstraintFieldGadget<E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(self))]
    fn to_constraint_field(&self) -> R1CSResult<Vec<FpVar<E2Fq<E>>>> {
        let mut res = Vec::new();

        let mut comm_gadget = self.comm.to_constraint_field()?;

        res.append(&mut comm_gadget);

        if self.shifted_comm.as_ref().is_some() {
            let mut shifted_comm_gadget =
                self.shifted_comm.as_ref().unwrap().to_constraint_field()?;
            res.append(&mut shifted_comm_gadget);
        }

        Ok(res)
    }
}

impl<E, PG> ToBytesGadget<E2Fq<E>> for CommitmentVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(self))]
    fn to_bytes(&self) -> R1CSResult<Vec<UInt8<E2Fq<E>>>> {
        let zero_shifted_comm = PG::G1Var::zero();

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.comm.to_bytes()?);

        let shifted_comm = self.shifted_comm.clone().unwrap_or(zero_shifted_comm);
        bytes.extend_from_slice(&shifted_comm.to_bytes()?);
        Ok(bytes)
    }
}

/// High level variable for a `MarlinKZG10` polynomial commitment, prepared for use in arirthmetic.
/// (`shifted_comm` is not prepared, due to the specific use case.)
#[derive(Derivative)]
#[derivative(Clone(bound = "E: PairingFriendlyCycle"))]
pub struct PreparedCommitmentVar<E: PairingFriendlyCycle, PG: PairingVar<E::Engine2, E2Fq<E>>> {
    prepared_comm: Vec<PG::G1Var>,
    shifted_comm: Option<PG::G1Var>,
}

impl<E, PG> PrepareGadget<CommitmentVar<E, PG>, E2Fq<E>> for PreparedCommitmentVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(unprepared))]
    fn prepare(unprepared: &CommitmentVar<E, PG>) -> R1CSResult<Self> {
        let mut prepared_comm = Vec::<PG::G1Var>::new();
        let supported_bits = <E2Fr<E> as PrimeField>::size_in_bits();

        let mut cur: PG::G1Var = unprepared.comm.clone();
        for _ in 0..supported_bits {
            prepared_comm.push(cur.clone());
            cur.double_in_place()?;
        }

        Ok(Self {
            prepared_comm,
            shifted_comm: unprepared.shifted_comm.clone(),
        })
    }
}

impl<E, PG> AllocVar<PreparedCommitment<E::Engine2>, E2Fq<E>> for PreparedCommitmentVar<E, PG>
where
    E: PairingFriendlyCycle,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable<T>(
        cs: impl Into<Namespace<E2Fq<E>>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> R1CSResult<Self>
    where
        T: Borrow<PreparedCommitment<E::Engine2>>,
    {
        let t = f()?;
        let obj = t.borrow();

        let ns = cs.into();
        let cs = ns.cs();

        let mut prepared_comm = Vec::<PG::G1Var>::new();

        for comm_elem in obj.prepared_comm.0.iter() {
            prepared_comm.push(<PG::G1Var as AllocVar<
                <E::Engine2 as PairingEngine>::G1Projective,
                E2Fq<E>,
            >>::new_variable(
                ark_relations::ns!(cs, "comm_elem"),
                || {
                    Ok(<<E::Engine2 as PairingEngine>::G1Projective as From<
                        <E::Engine2 as PairingEngine>::G1Affine,
                    >>::from(*comm_elem))
                },
                mode,
            )?);
        }

        let shifted_comm = if obj.shifted_comm.is_some() {
            Some(<PG::G1Var as AllocVar<
                <E::Engine2 as PairingEngine>::G1Projective,
                E2Fq<E>,
            >>::new_variable(
                ark_relations::ns!(cs, "shifted_comm"),
                || {
                    Ok(<<E::Engine2 as PairingEngine>::G1Projective as From<
                        <E::Engine2 as PairingEngine>::G1Affine,
                    >>::from(obj.shifted_comm.unwrap().0))
                },
                mode,
            )?)
        } else {
            None
        };

        Ok(Self {
            prepared_comm,
            shifted_comm,
        })
    }
}

/// High level variable for a `MarlinKZG10` polynomial commitment, along with a string label and a
/// degree bound.
#[derive(Derivative)]
#[derivative(Clone(bound = "E: PairingFriendlyCycle"))]
pub struct LabeledCommitmentVar<E: PairingFriendlyCycle, PG: PairingVar<E::Engine2, E2Fq<E>>>
where
    E2Fq<E>: PrimeField,
{
    /// A text label for the commitment.
    pub label: String,
    /// The plain commitment.
    pub commitment: CommitmentVar<E, PG>,
    /// Optionally, a bound on the polynomial degree.
    pub degree_bound: Option<FpVar<E2Fq<E>>>,
}

impl<E, PG> AllocVar<LabeledCommitment<Commitment<E::Engine2>>, E2Fq<E>>
    for LabeledCommitmentVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, value_gen))]
    fn new_variable<T>(
        cs: impl Into<Namespace<E2Fq<E>>>,
        value_gen: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> R1CSResult<Self>
    where
        T: Borrow<LabeledCommitment<Commitment<E::Engine2>>>,
    {
        value_gen().and_then(|labeled_commitment| {
            let ns = cs.into();
            let cs = ns.cs();

            let labeled_commitment = labeled_commitment.borrow().clone();
            let label = labeled_commitment.label().to_string();
            let commitment = labeled_commitment.commitment();
            let degree_bound = labeled_commitment.degree_bound();

            let commitment = CommitmentVar::new_variable(
                ark_relations::ns!(cs, "commitment"),
                || Ok(commitment),
                mode,
            )?;

            let degree_bound = if let Some(degree_bound) = degree_bound {
                FpVar::<E2Fq<E>>::new_variable(
                    ark_relations::ns!(cs, "degree_bound"),
                    || Ok(<E2Fq<E> as From<u128>>::from(degree_bound as u128)),
                    mode,
                )
                .ok()
            } else {
                None
            };

            Ok(Self {
                label,
                commitment,
                degree_bound,
            })
        })
    }
}

/// High level variable for a `MarlinKZG10` polynomial commitment, along with a string label and a
/// degree bound, prepared for use in arithmetic.
#[derive(Derivative)]
#[derivative(Clone(bound = "E: PairingFriendlyCycle"))]
pub struct PreparedLabeledCommitmentVar<
    E: PairingFriendlyCycle,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
> where
    E2Fq<E>: PrimeField,
{
    /// A text label for the commitment.
    pub label: String,
    /// The plain commitment.
    pub prepared_commitment: PreparedCommitmentVar<E, PG>,
    /// Optionally, a bound on the polynomial degree.
    pub degree_bound: Option<FpVar<E2Fq<E>>>,
}

impl<E, PG> PrepareGadget<LabeledCommitmentVar<E, PG>, E2Fq<E>>
    for PreparedLabeledCommitmentVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(unprepared))]
    fn prepare(unprepared: &LabeledCommitmentVar<E, PG>) -> R1CSResult<Self> {
        let prepared_commitment = PreparedCommitmentVar::prepare(&unprepared.commitment)?;

        Ok(Self {
            label: unprepared.label.clone(),
            prepared_commitment,
            degree_bound: unprepared.degree_bound.clone(),
        })
    }
}

/// High level variable for a `MarlinKZG10` opening proof.
#[allow(clippy::type_complexity)]
#[derive(Derivative)]
#[derivative(Clone(bound = "E: PairingFriendlyCycle"))]
pub struct ProofVar<E: PairingFriendlyCycle, PG: PairingVar<E::Engine2, E2Fq<E>>>
where
    E2Fq<E>: PrimeField,
{
    /// The commitment to the witness polynomial.
    pub w: PG::G1Var,
    /// The evaluation of the random hiding polynomial.
    pub random_v: Option<NonNativeFieldVar<E2Fr<E>, E2Fq<E>>>,
}

impl<E, PG> AllocVar<Proof<E::Engine2>, E2Fq<E>> for ProofVar<E, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, value_gen))]
    fn new_variable<T>(
        cs: impl Into<Namespace<E2Fq<E>>>,
        value_gen: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> R1CSResult<Self>
    where
        T: Borrow<Proof<E::Engine2>>,
    {
        value_gen().and_then(|proof| {
            let ns = cs.into();
            let cs = ns.cs();

            let Proof { w, random_v } = *proof.borrow();
            let w = PG::G1Var::new_variable(ark_relations::ns!(cs, "w"), || Ok(w), mode)?;

            let random_v = match random_v {
                None => None,
                Some(random_v_inner) => Some(NonNativeFieldVar::new_variable(
                    ark_relations::ns!(cs, "random_v"),
                    || Ok(random_v_inner),
                    mode,
                )?),
            };

            Ok(Self { w, random_v })
        })
    }
}

/// High level variable for a batched `MarlinKZG10` proof, for the opening of a linear combination
/// of polynomials.
#[allow(clippy::type_complexity)]
#[derive(Derivative)]
#[derivative(Clone(bound = "E: PairingFriendlyCycle"))]
pub struct BatchLCProofVar<
    E: PairingFriendlyCycle,
    P: UVPolynomial<E2Fr<E>, Point = E2Fr<E>>,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
> where
    E2Fq<E>: PrimeField,
{
    /// Evaluation proofs.
    pub proofs: Vec<ProofVar<E, PG>>,
    /// Evaluations required to verify the proof.
    pub evals: Option<Vec<NonNativeFieldVar<E2Fr<E>, E2Fq<E>>>>,
    #[doc(hidden)]
    pub polynomial: PhantomData<P>,
}

impl<E, P, PG> AllocVar<BatchLCProof<E2Fr<E>, P, MarlinKZG10<E::Engine2, P>>, E2Fq<E>>
    for BatchLCProofVar<E, P, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    P: UVPolynomial<E2Fr<E>, Point = E2Fr<E>>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, value_gen))]
    fn new_variable<T>(
        cs: impl Into<Namespace<E2Fq<E>>>,
        value_gen: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> R1CSResult<Self>
    where
        T: Borrow<BatchLCProof<E2Fr<E>, P, MarlinKZG10<E::Engine2, P>>>,
    {
        value_gen().map(|proof| {
            let ns = cs.into();
            let cs = ns.cs();

            let BatchLCProof { proof, evals } = proof.borrow().clone();

            let proofs: Vec<Proof<_>> = proof.to_vec();
            let proofs: Vec<ProofVar<E, PG>> = proofs
                .iter()
                .map(|p| {
                    ProofVar::new_variable(ark_relations::ns!(cs, "proof"), || Ok(p), mode).unwrap()
                })
                .collect();

            #[allow(clippy::type_complexity)]
            let evals: Option<Vec<NonNativeFieldVar<E2Fr<E>, E2Fq<E>>>> = match evals {
                None => None,
                Some(evals_inner) => Some(
                    evals_inner
                        .iter()
                        .map(|e| {
                            NonNativeFieldVar::new_variable(
                                ark_relations::ns!(cs, "evaluation"),
                                || Ok(e),
                                mode,
                            )
                            .unwrap()
                        })
                        .collect(),
                ),
            };

            Self {
                proofs,
                evals,
                polynomial: PhantomData,
            }
        })
    }
}

/// Helper struct for information about one member of a linear combination.
#[derive(Derivative)]
#[derivative(Clone(bound = "E: PairingFriendlyCycle"))]
pub struct LCItem<E, PG>
where
    E: PairingFriendlyCycle,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
    E2Fq<E>: PrimeField,
{
    coeff: Option<NonNativeFieldVar<E2Fr<E>, E2Fq<E>>>,
    degree_bound: Option<FpVar<E2Fq<E>>>,
    comm: PreparedCommitmentVar<E, PG>,
    negate: bool,
}

/// Gadget for the `MarlinKZG10` polynomial commitment verifier.
#[derive(Derivative)]
#[derivative(Clone(bound = "E: PairingFriendlyCycle"))]
pub struct MarlinKZG10Gadget<E, P, PG>
where
    E: PairingFriendlyCycle,
    P: UVPolynomial<E2Fr<E>, Point = E2Fr<E>>,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    _cycle_engine: PhantomData<E>,
    _pairing_gadget: PhantomData<PG>,
    _polynomial: PhantomData<P>,
}

impl<E, P, PG> MarlinKZG10Gadget<E, P, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    P: UVPolynomial<E2Fr<E>, Point = E2Fr<E>>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    #[allow(clippy::type_complexity, clippy::too_many_arguments)]
    #[tracing::instrument(
        target = "r1cs",
        skip(prepared_verification_key, lc_info, query_set, evaluations, proofs)
    )]
    fn prepared_batch_check_evaluations(
        cs: ConstraintSystemRef<E2Fq<E>>,
        prepared_verification_key: &<Self as PCCheckVar<
            E2Fr<E>,
            P,
            MarlinKZG10<E::Engine2, P>,
            E2Fq<E>,
        >>::PreparedVerifierKeyVar,
        lc_info: &[(
            String,
            Vec<LCItem<E, PG>>,
        )],
        query_set: &QuerySetVar<E2Fr<E>, E2Fq<E>>,
        evaluations: &EvaluationsVar<E2Fr<E>, E2Fq<E>>,
        proofs: &[<Self as PCCheckVar<
            E2Fr<E>,
            P,
            MarlinKZG10<E::Engine2, P>,
            E2Fq<E>,
        >>::ProofVar],
        opening_challenges: &[NonNativeFieldVar<E2Fr<E>, E2Fq<E>>],
        opening_challenges_bits: &[Vec<Boolean<E2Fq<E>>>],
        batching_rands: &[NonNativeFieldVar<E2Fr<E>, E2Fq<E>>],
        batching_rands_bits: &[Vec<Boolean<E2Fq<E>>>],
    ) -> R1CSResult<Boolean<E2Fq<E>>> {
        let mut batching_rands = batching_rands.to_vec();
        let mut batching_rands_bits = batching_rands_bits.to_vec();

        let commitment_lcs: BTreeMap<
            String,
            (
                String,
                Vec<LCItem<E, PG>>,
            ),
        > = lc_info.iter().map(|c| (c.0.clone(), c.clone())).collect();

        let mut query_to_labels_map = BTreeMap::new();

        for (label, point) in query_set.0.iter() {
            let labels = query_to_labels_map
                .entry(point.name.clone())
                .or_insert((point.value.clone(), BTreeSet::new()));
            labels.1.insert(label);
        }

        // Accumulate commitments and evaluations for each query.
        let mut combined_queries = Vec::new();
        let mut combined_comms = Vec::new();
        let mut combined_evals = Vec::new();
        for (_, (point, labels)) in query_to_labels_map.into_iter() {
            let mut comms_to_combine = Vec::<
                Vec<LCItem<E, PG>>,
            >::new();
            let mut values_to_combine = Vec::new();
            for label in labels.into_iter() {
                let commitment_lc = commitment_lcs.get(label).unwrap().clone();

                let v_i = evaluations
                    .0
                    .get(&LabeledPointVar {
                        name: label.clone(),
                        value: point.clone(),
                    })
                    .unwrap();

                comms_to_combine.push(commitment_lc.1.clone());
                values_to_combine.push(v_i.clone());
            }

            // Accumulate the commitments and evaluations corresponding to `query`.
            let mut combined_comm = PG::G1Var::zero();
            let mut combined_eval = NonNativeFieldMulResultVar::<E2Fr<E>, E2Fq<E>>::zero();

            let mut opening_challenges_counter = 0;

            for (commitment_lcs, value) in comms_to_combine.into_iter().zip(values_to_combine) {
                let challenge = opening_challenges[opening_challenges_counter].clone();

                let challenge_bits = opening_challenges_bits[opening_challenges_counter].clone();
                opening_challenges_counter += 1;

                for lc_info in commitment_lcs.iter() {
                    let LCItem::<E, PG> { coeff, degree_bound, comm, negate } = lc_info;
                    let PreparedCommitmentVar { shifted_comm, .. } = comm;

                    if coeff.is_none() {
                        // To combine the commitments, we multiply each by one of the random challenges, and sum.
                        let mut comm_times_challenge = PG::G1Var::zero();
                        {
                            for (bit, base_power) in
                                challenge_bits.iter().zip(comm.prepared_comm.iter())
                            {
                                let mut new_encoded = base_power.clone();
                                new_encoded += comm_times_challenge.clone();
                                comm_times_challenge = PG::G1Var::conditionally_select(
                                    bit,
                                    &new_encoded,
                                    &comm_times_challenge,
                                )?;
                            }
                        }

                        if negate.eq(&true) {
                            comm_times_challenge = comm_times_challenge.negate()?;
                        }

                        combined_comm += comm_times_challenge.clone();

                        // If the degree bound is specified, we include the adjusted degree-shifted commitment
                        // (that is, c_i' - v_i beta^{D - d_i} G), where d_i is the specific degree bound and
                        // v_i is the evaluation, in the combined commitment,
                        if let Some(degree_bound) = degree_bound {
                            let challenge_shifted_bits =
                                opening_challenges_bits[opening_challenges_counter].clone();
                            opening_challenges_counter += 1;

                            let mut shifted_comm = shifted_comm.clone().unwrap();

                            if negate.eq(&true) {
                                shifted_comm = shifted_comm.negate()?;
                            }

                            let value_bits = value.to_bits_le()?;
                            let shift_power = prepared_verification_key
                                .get_shift_power(cs.clone(), degree_bound)
                                .unwrap();

                            let mut shift_power_times_value = PG::G1Var::zero();
                            {
                                for (bit, base_power) in value_bits.iter().zip(&shift_power) {
                                    let mut new_encoded = base_power.clone();
                                    new_encoded += shift_power_times_value.clone();
                                    shift_power_times_value = PG::G1Var::conditionally_select(
                                        bit,
                                        &new_encoded,
                                        &shift_power_times_value,
                                    )?;
                                }
                            }
                            let mut adjusted_comm = shifted_comm;
                            adjusted_comm -= shift_power_times_value;
                            let adjusted_comm_times_challenge =
                                adjusted_comm.scalar_mul_le(challenge_shifted_bits.iter())?;
                            combined_comm += adjusted_comm_times_challenge;
                        }
                    } else {
                        assert!(degree_bound.is_none());

                        let mut comm_times_challenge = PG::G1Var::zero();
                        let coeff = coeff.clone().unwrap();

                        let challenge_times_coeff = &challenge * &coeff;
                        let challenge_times_coeff_bits = challenge_times_coeff.to_bits_le()?;

                        {
                            for (bit, base_power) in
                                challenge_times_coeff_bits.iter().zip(&comm.prepared_comm)
                            {
                                let mut new_encoded = comm_times_challenge.clone();
                                new_encoded += base_power.clone();
                                comm_times_challenge = PG::G1Var::conditionally_select(
                                    bit,
                                    &new_encoded,
                                    &comm_times_challenge,
                                )?;
                            }
                        }

                        if negate.eq(&true) {
                            comm_times_challenge = comm_times_challenge.negate()?;
                        }

                        combined_comm += comm_times_challenge;
                    }
                }
                // Similarly, we add up the evaluations, multiplied with random challenges.
                let value_times_challenge_unreduced = value.mul_without_reduce(&challenge)?;

                combined_eval += &value_times_challenge_unreduced;
            }

            let combined_eval_reduced = combined_eval.reduce()?;

            combined_queries.push(point.clone());
            combined_comms.push(combined_comm);
            combined_evals.push(combined_eval_reduced);
        }

        // Perform the batch check.
        {
            let mut total_c = PG::G1Var::zero();
            let mut total_w = PG::G1Var::zero();

            let mut g_multiplier = NonNativeFieldMulResultVar::<E2Fr<E>, E2Fq<E>>::zero();
            let mut g_multiplier_reduced = NonNativeFieldVar::<E2Fr<E>, E2Fq<E>>::zero();
            for (i, (((c, z), v), proof)) in combined_comms
                .iter()
                .zip(combined_queries)
                .zip(combined_evals)
                .zip(proofs)
                .enumerate()
            {
                let z_bits = z.to_bits_le()?;

                let w_times_z = proof.w.scalar_mul_le(z_bits.iter())?;

                let mut c_plus_w_times_z = c.clone();
                c_plus_w_times_z += w_times_z;

                if i != 0 {
                    let randomizer = batching_rands.remove(0);
                    let randomizer_bits = batching_rands_bits.remove(0);

                    let randomizer_times_v = randomizer.mul_without_reduce(&v)?;

                    g_multiplier += &randomizer_times_v;

                    let c_times_randomizer =
                        c_plus_w_times_z.scalar_mul_le(randomizer_bits.iter())?;
                    let w_times_randomizer = proof.w.scalar_mul_le(randomizer_bits.iter())?;
                    total_c += c_times_randomizer;
                    total_w += w_times_randomizer;
                } else {
                    g_multiplier_reduced += v;
                    total_c += c_plus_w_times_z;
                    total_w += proof.w.clone();
                }
            }

            // Prepare each input to the pairing.
            let (prepared_total_w, prepared_beta_h, prepared_total_c, prepared_h) = {
                let g_multiplier_reduced = g_multiplier.reduce()? + &g_multiplier_reduced;
                let g_multiplier_bits = g_multiplier_reduced.to_bits_le()?;

                let mut g_times_mul = PG::G1Var::zero();
                {
                    for (bit, base_power) in g_multiplier_bits
                        .iter()
                        .zip(&prepared_verification_key.prepared_g)
                    {
                        let mut new_encoded = g_times_mul.clone();
                        new_encoded += base_power.clone();
                        g_times_mul =
                            PG::G1Var::conditionally_select(bit, &new_encoded, &g_times_mul)?;
                    }
                }
                total_c -= g_times_mul;
                total_w = total_w.negate()?;

                let prepared_total_w = PG::prepare_g1(&total_w)?;
                let prepared_beta_h = prepared_verification_key.prepared_beta_h.clone();
                let prepared_total_c = PG::prepare_g1(&total_c)?;
                let prepared_h = prepared_verification_key.prepared_h.clone();

                (
                    prepared_total_w,
                    prepared_beta_h,
                    prepared_total_c,
                    prepared_h,
                )
            };

            let lhs = PG::product_of_pairings(
                &[prepared_total_w, prepared_total_c],
                &[prepared_beta_h, prepared_h],
            )?;

            let rhs = &PG::GTVar::one();
            lhs.is_eq(&rhs)
        }
    }
}

impl<E, P, PG> PCCheckVar<E2Fr<E>, P, MarlinKZG10<E::Engine2, P>, E2Fq<E>>
    for MarlinKZG10Gadget<E, P, PG>
where
    E: PairingFriendlyCycle,
    E2Fq<E>: PrimeField,
    P: UVPolynomial<E2Fr<E>, Point = E2Fr<E>>,
    for<'a, 'b> &'a P: Div<&'b P, Output = P>,
    PG: PairingVar<E::Engine2, E2Fq<E>>,
{
    type VerifierKeyVar = VerifierKeyVar<E, PG>;
    type PreparedVerifierKeyVar = PreparedVerifierKeyVar<E, PG>;
    type CommitmentVar = CommitmentVar<E, PG>;
    type PreparedCommitmentVar = PreparedCommitmentVar<E, PG>;
    type LabeledCommitmentVar = LabeledCommitmentVar<E, PG>;
    type PreparedLabeledCommitmentVar = PreparedLabeledCommitmentVar<E, PG>;
    type ProofVar = ProofVar<E, PG>;
    type BatchLCProofVar = BatchLCProofVar<E, P, PG>;

    #[allow(clippy::type_complexity)]
    #[tracing::instrument(
        target = "r1cs",
        skip(verification_key, commitments, query_set, evaluations, proofs)
    )]
    fn batch_check_evaluations(
        _cs: ConstraintSystemRef<E2Fq<E>>,
        verification_key: &Self::VerifierKeyVar,
        commitments: &[Self::LabeledCommitmentVar],
        query_set: &QuerySetVar<E2Fr<E>, E2Fq<E>>,
        evaluations: &EvaluationsVar<E2Fr<E>, E2Fq<E>>,
        proofs: &[Self::ProofVar],
        rand_data: &PCCheckRandomDataVar<E2Fr<E>, E2Fq<E>>,
    ) -> R1CSResult<Boolean<E2Fq<E>>> {
        let mut batching_rands = rand_data.batching_rands.to_vec();
        let mut batching_rands_bits = rand_data.batching_rands_bits.to_vec();

        let commitments: BTreeMap<_, _> =
            commitments.iter().map(|c| (c.label.clone(), c)).collect();
        let mut query_to_labels_map = BTreeMap::new();

        for (label, point) in query_set.0.iter() {
            let labels = query_to_labels_map
                .entry(point.name.clone())
                .or_insert((point.value.clone(), BTreeSet::new()));
            labels.1.insert(label);
        }

        // Accumulate commitments and evaluations for each query.
        let mut combined_queries = Vec::new();
        let mut combined_comms = Vec::new();
        let mut combined_evals = Vec::new();
        for (_, (point, labels)) in query_to_labels_map.into_iter() {
            let mut comms_to_combine: Vec<Self::LabeledCommitmentVar> = Vec::new();
            let mut values_to_combine = Vec::new();
            for label in labels.into_iter() {
                let commitment = &(*commitments.get(label).unwrap()).clone();
                let degree_bound = commitment.degree_bound.clone();
                assert_eq!(
                    degree_bound.is_some(),
                    commitment.commitment.shifted_comm.is_some()
                );

                let v_i = evaluations
                    .0
                    .get(&LabeledPointVar {
                        name: label.clone(),
                        value: point.clone(),
                    })
                    .unwrap();

                comms_to_combine.push(commitment.clone());
                values_to_combine.push(v_i.clone());
            }

            // Accumulate the commitments and evaluations corresponding to `query`.
            let mut combined_comm = PG::G1Var::zero();
            let mut combined_eval = NonNativeFieldMulResultVar::<E2Fr<E>, E2Fq<E>>::zero();

            let mut opening_challenges_counter = 0;

            for (labeled_commitment, value) in
                comms_to_combine.into_iter().zip(values_to_combine.iter())
            {
                let challenge = rand_data.opening_challenges[opening_challenges_counter].clone();
                let challenge_bits =
                    rand_data.opening_challenges_bits[opening_challenges_counter].clone();
                opening_challenges_counter += 1;

                let LabeledCommitmentVar {
                    commitment,
                    degree_bound,
                    ..
                } = labeled_commitment;
                let CommitmentVar { shifted_comm, .. } = commitment;

                // To combine the commitments, we multiply each by one of the random challenges, and sum.
                combined_comm += commitment.comm.scalar_mul_le(challenge_bits.iter())?;

                // Similarly, we add up the evaluations, multiplied with random challenges.
                let value_times_challenge_unreduced = value.mul_without_reduce(&challenge)?;
                combined_eval += &value_times_challenge_unreduced;

                // If the degree bound is specified, we include the adjusted degree-shifted commitment
                // (that is, c_i' - v_i beta^{D - d_i} G), where d_i is the specific degree bound and
                // v_i is the evaluation, in the cocmbined commitment,
                if let Some(degree_bound) = degree_bound {
                    let challenge_shifted_bits =
                        rand_data.opening_challenges_bits[opening_challenges_counter].clone();
                    opening_challenges_counter += 1;

                    let shifted_comm = shifted_comm.as_ref().unwrap().clone();

                    let value_bits = value.to_bits_le()?;
                    let shift_power = verification_key
                        .get_shift_power(verification_key.g.cs(), &degree_bound)
                        .unwrap();

                    let shift_power_times_value = shift_power.scalar_mul_le(value_bits.iter())?;
                    let mut adjusted_comm = shifted_comm;
                    adjusted_comm -= shift_power_times_value;

                    let adjusted_comm_times_challenge =
                        adjusted_comm.scalar_mul_le(challenge_shifted_bits.iter())?;

                    combined_comm += adjusted_comm_times_challenge;
                }
            }

            combined_queries.push(point.clone());
            combined_comms.push(combined_comm);
            combined_evals.push(combined_eval);
        }

        // Perform the batch check.
        {
            let mut total_c = PG::G1Var::zero();
            let mut total_w = PG::G1Var::zero();

            let mut g_multiplier = NonNativeFieldMulResultVar::<E2Fr<E>, E2Fq<E>>::zero();
            for (((c, z), v), proof) in combined_comms
                .iter()
                .zip(combined_queries)
                .zip(combined_evals)
                .zip(proofs)
            {
                let z_bits = z.to_bits_le()?;

                let w_times_z = proof.w.scalar_mul_le(z_bits.iter())?;
                let mut c_plus_w_times_z = c.clone();
                c_plus_w_times_z += w_times_z;

                let randomizer = batching_rands.remove(0);
                let randomizer_bits = batching_rands_bits.remove(0);

                let v_reduced = v.reduce()?;
                let randomizer_times_v = randomizer.mul_without_reduce(&v_reduced)?;

                g_multiplier += randomizer_times_v;

                let c_times_randomizer = c_plus_w_times_z.scalar_mul_le(randomizer_bits.iter())?;
                let w_times_randomizer = proof.w.scalar_mul_le(randomizer_bits.iter())?;
                total_c += c_times_randomizer;
                total_w += w_times_randomizer;
            }

            // Prepare each input to the pairing.
            let (prepared_total_w, prepared_beta_h, prepared_total_c, prepared_h) = {
                let g_multiplier_reduced = g_multiplier.reduce()?;
                let g_multiplier_bits = g_multiplier_reduced.to_bits_le()?;

                let g_times_mul = verification_key.g.scalar_mul_le(g_multiplier_bits.iter())?;
                total_c -= g_times_mul;
                total_w = total_w.negate()?;

                let prepared_total_w = PG::prepare_g1(&total_w)?;
                let prepared_beta_h = PG::prepare_g2(&verification_key.beta_h)?;
                let prepared_total_c = PG::prepare_g1(&total_c)?;
                let prepared_h = PG::prepare_g2(&verification_key.h)?;

                (
                    prepared_total_w,
                    prepared_beta_h,
                    prepared_total_c,
                    prepared_h,
                )
            };

            let lhs = PG::product_of_pairings(
                &[prepared_total_w, prepared_total_c],
                &[prepared_beta_h, prepared_h],
            )?;

            let rhs = &PG::GTVar::one();

            lhs.is_eq(rhs)
        }
    }

    #[allow(clippy::type_complexity)]
    #[tracing::instrument(
        target = "r1cs",
        skip(
            prepared_verification_key,
            linear_combinations,
            prepared_commitments,
            query_set,
            proof,
            evaluations
        )
    )]
    fn prepared_check_combinations(
        cs: ConstraintSystemRef<E2Fq<E>>,
        prepared_verification_key: &Self::PreparedVerifierKeyVar,
        linear_combinations: &[LinearCombinationVar<E2Fr<E>, E2Fq<E>>],
        prepared_commitments: &[Self::PreparedLabeledCommitmentVar],
        query_set: &QuerySetVar<E2Fr<E>, E2Fq<E>>,
        evaluations: &EvaluationsVar<E2Fr<E>, E2Fq<E>>,
        proof: &Self::BatchLCProofVar,
        rand_data: &PCCheckRandomDataVar<E2Fr<E>, E2Fq<E>>,
    ) -> R1CSResult<Boolean<E2Fq<E>>> {
        let BatchLCProofVar { proofs, .. } = proof;

        let label_comm_map = prepared_commitments
            .iter()
            .map(|c| (c.label.clone(), c))
            .collect::<BTreeMap<_, _>>();

        let mut lc_info = Vec::new();
        let mut evaluations = evaluations.clone();

        // For each linear combination, we sum up the relevant commitments, multiplied
        // with their corresponding coefficients; these combined commitments are then
        // the inputs to the normal batch check.
        for lc in linear_combinations.iter() {
            let lc_label = lc.label.clone();
            let num_polys = lc.terms.len();

            let mut coeffs_and_comms = Vec::new();

            for (coeff, label) in lc.terms.iter() {
                if label.is_one() {
                    for (label, ref mut eval) in evaluations.0.iter_mut() {
                        if label.name == lc_label {
                            match coeff.clone() {
                                LinearCombinationCoeffVar::One => {
                                    **eval = (**eval).clone() - &NonNativeFieldVar::one()
                                }
                                LinearCombinationCoeffVar::MinusOne => {
                                    **eval = (**eval).clone() + &NonNativeFieldVar::one()
                                }
                                LinearCombinationCoeffVar::Var(variable) => {
                                    **eval = (**eval).clone() - &variable
                                }
                            };
                        }
                    }
                } else {
                    let label: &String = label.try_into().unwrap();
                    let &cur_comm = label_comm_map.get(label).unwrap();
                    let negate = match coeff {
                        LinearCombinationCoeffVar::One | LinearCombinationCoeffVar::Var(_) => false,
                        LinearCombinationCoeffVar::MinusOne => true,
                    };

                    if num_polys == 1 && cur_comm.degree_bound.is_some() {
                        assert!(!negate);
                    }

                    let coeff = match coeff {
                        LinearCombinationCoeffVar::One => None,
                        LinearCombinationCoeffVar::MinusOne => None,
                        LinearCombinationCoeffVar::Var(variable) => Some(variable.clone()),
                    };

                    coeffs_and_comms.push( LCItem {
                        coeff,
                        degree_bound: cur_comm.degree_bound.clone(),
                        comm: cur_comm.prepared_commitment.clone(),
                        negate
                    });
                }
            }

            lc_info.push((lc_label, coeffs_and_comms));
        }

        Self::prepared_batch_check_evaluations(
            cs,
            prepared_verification_key,
            lc_info.as_slice(),
            &query_set,
            &evaluations,
            proofs,
            &rand_data.opening_challenges,
            &rand_data.opening_challenges_bits,
            &rand_data.batching_rands,
            &rand_data.batching_rands_bits,
        )
    }

    fn create_labeled_commitment(
        label: String,
        commitment: Self::CommitmentVar,
        degree_bound: Option<FpVar<E2Fq<E>>>,
    ) -> Self::LabeledCommitmentVar {
        Self::LabeledCommitmentVar {
            label,
            commitment,
            degree_bound,
        }
    }

    fn create_prepared_labeled_commitment(
        label: String,
        prepared_commitment: Self::PreparedCommitmentVar,
        degree_bound: Option<FpVar<E2Fq<E>>>,
    ) -> Self::PreparedLabeledCommitmentVar {
        Self::PreparedLabeledCommitmentVar {
            label,
            prepared_commitment,
            degree_bound,
        }
    }
}
