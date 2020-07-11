use crate::{
    data_structures::{LCTerm, LabeledCommitment, LinearCombination},
    kzg10::{Proof, VerifierKey as KZG10VerifierKey},
    marlin_pc::{
        data_structures::{Commitment, VerifierKey},
        MarlinKZG10, PreparedCommitment, PreparedVerifierKey,
    },
    pc_constraints::{EvaluationsGadget, PCCheckGadget, QuerySetGadget},
    BTreeMap,
    BTreeSet,
    BatchLCProof, // PolynomialCommitment,
    PrepareGadget,
    String,
    ToString,
    Vec,
};
use algebra::{AffineCurve, CycleEngine, PairingEngine};
use algebra_core::{PrimeField, ProjectiveCurve};
use core::{borrow::Borrow, convert::TryInto, marker::PhantomData, ops::MulAssign};
use nonnative::{NonNativeFieldGadget, NonNativeFieldMulResultGadget, NonNativeFieldParams};
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::{
    alloc::AllocGadget,
    bits::{boolean::Boolean, uint32::UInt32, uint8::UInt8},
    eq::EqGadget,
    fields::{fp::FpGadget, FieldGadget, ToConstraintFieldGadget},
    groups::GroupGadget,
    pairing::PairingGadget,
    select::CondSelectGadget,
    ToBitsGadget, ToBytesGadget,
};

/// Gadget for the verification key of the Marlin-KZG10 polynomial commitment scheme.
pub struct VerifierKeyGadget<
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
> where
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    /// Generator of G1.
    pub g: PG::G1Gadget,
    /// Generator of G2.
    pub h: PG::G2Gadget,
    /// Generator of G1, times first monomial.
    pub beta_h: PG::G2Gadget,
    /// Used for the shift powers associated with different degree bounds.
    pub degree_bounds_and_shift_powers: Option<Vec<(usize, PG::G1Gadget)>>,
    np_phantom: PhantomData<NP>,
}

impl<CycleE, NP, PG> VerifierKeyGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    /// Find the appropriate shift for the degree bound.
    pub fn get_shift_power(&self, bound: usize) -> Option<PG::G1Gadget> {
        self.degree_bounds_and_shift_powers.as_ref().and_then(|v| {
            v.binary_search_by(|(d, _)| d.cmp(&bound))
                .ok()
                .map(|i| v[i].1.clone())
        })
    }
}

impl<CycleE, NP, PG> Clone for VerifierKeyGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn clone(&self) -> Self {
        VerifierKeyGadget {
            g: self.g.clone(),
            h: self.h.clone(),
            beta_h: self.beta_h.clone(),
            degree_bounds_and_shift_powers: self.degree_bounds_and_shift_powers.clone(),
            np_phantom: PhantomData,
        }
    }
}

impl<CycleE, NP, PG> AllocGadget<VerifierKey<CycleE::E2>, <CycleE::E1 as PairingEngine>::Fr>
    for VerifierKeyGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn alloc_constant<T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        _cs: CS,
        _val: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<VerifierKey<CycleE::E2>>,
    {
        unimplemented!()
    }

    #[inline]
    fn alloc<FN, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<VerifierKey<CycleE::E2>>,
    {
        value_gen().and_then(|vk| {
            let vk_orig = vk.borrow().clone();

            let VerifierKey {
                vk,
                degree_bounds_and_shift_powers,
                ..
            } = vk_orig;

            let degree_bounds_and_shift_powers = degree_bounds_and_shift_powers
                .and_then(|vec| Some(vec.iter()
                                        .map(|(s, g)|
                                        //(UInt32::alloc(cs.ns(|| format!("map {}", s)), Some(*s as u32)),
                                        (*s,
                                        PG::G1Gadget::alloc(cs.ns(|| format!("pow for {}", s)), || Ok(g.into_projective())).unwrap())
                                ).collect()));

            let KZG10VerifierKey { g, h, beta_h, .. } = vk;

            let g = PG::G1Gadget::alloc(cs.ns(|| "g"), || Ok(g.into_projective()))?;
            let h = PG::G2Gadget::alloc(cs.ns(|| "h"), || Ok(h.into_projective()))?;
            let beta_h = PG::G2Gadget::alloc(cs.ns(|| "beta_h"), || Ok(beta_h.into_projective()))?;

            Ok(Self { g, h, beta_h, degree_bounds_and_shift_powers, np_phantom: PhantomData })
        })
    }

    #[inline]
    fn alloc_input<FN, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<VerifierKey<CycleE::E2>>,
    {
        value_gen().and_then(|vk| {
            let vk_orig = vk.borrow().clone();

            let VerifierKey {
                vk,
                degree_bounds_and_shift_powers,
                ..
            } = vk_orig;

            let degree_bounds_and_shift_powers = degree_bounds_and_shift_powers
                .and_then(|vec| Some(vec.iter()
                                        .map(|(s, g)|
                                        //(UInt32::alloc(cs.ns(|| format!("map {}", s)), Some(*s as u32)),
                                        (*s,
                                        PG::G1Gadget::alloc(cs.ns(|| format!("pow for {}", s)), || Ok(g.into_projective())).unwrap())
                                ).collect()));

            let KZG10VerifierKey { g, h, beta_h, .. } = vk;
            let g = PG::G1Gadget::alloc(cs.ns(|| "g"), || Ok(g.into_projective()))?;
            let h = PG::G2Gadget::alloc(cs.ns(|| "h"), || Ok(h.into_projective()))?;
            let beta_h = PG::G2Gadget::alloc(cs.ns(|| "beta_h"), || Ok(beta_h.into_projective()))?;

            Ok(Self { g, h, beta_h, degree_bounds_and_shift_powers, np_phantom: PhantomData })
        })
    }
}

impl<CycleE, NP, PG> ToBytesGadget<<CycleE::E1 as PairingEngine>::Fr>
    for VerifierKeyGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    #[inline]
    fn to_bytes<CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.g.to_bytes(&mut cs.ns(|| "g to bytes"))?);
        bytes.extend_from_slice(&self.h.to_bytes(&mut cs.ns(|| "h to bytes"))?);
        bytes.extend_from_slice(&self.beta_h.to_bytes(&mut cs.ns(|| "beta_h to bytes"))?);

        Ok(bytes)
    }
}

impl<CycleE, NP, PG> ToConstraintFieldGadget<<CycleE::E1 as PairingEngine>::Fr>
    for VerifierKeyGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn to_field_gadgets<CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<FpGadget<<CycleE::E1 as PairingEngine>::Fr>>, SynthesisError> {
        let mut res = Vec::new();

        let mut g_gadget = self.g.to_field_gadgets(&mut cs.ns(|| "g"))?;
        let mut h_gadget = self.h.to_field_gadgets(&mut cs.ns(|| "h"))?;
        let mut beta_h_gadget = self.beta_h.to_field_gadgets(&mut cs.ns(|| "beta_h"))?;

        res.append(&mut g_gadget);
        res.append(&mut h_gadget);
        res.append(&mut beta_h_gadget);

        if self.degree_bounds_and_shift_powers.as_ref().is_some() {
            let list = self.degree_bounds_and_shift_powers.as_ref().unwrap();
            for (i, (d, shift_power)) in list.iter().enumerate() {
                let d_bigint = (d.clone() as u64).into();
                let d_gadget = FpGadget::<<CycleE::E1 as PairingEngine>::Fr>::alloc_constant(
                    &mut cs.ns(|| format!("d#{}", i)),
                    <CycleE::E1 as PairingEngine>::Fr::from_repr(d_bigint).unwrap(),
                )?;
                let mut shift_power_gadget =
                    shift_power.to_field_gadgets(&mut cs.ns(|| format!("shift_power#{}", i)))?;

                res.push(d_gadget);
                res.append(&mut shift_power_gadget);
            }
        }

        Ok(res)
    }
}

/// Gadget for the verification key of the Marlin-KZG10 polynomial commitment scheme.
pub struct PreparedVerifierKeyGadget<
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
> where
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    /// Generator of G1.
    pub prepared_g: Vec<PG::G1Gadget>,
    /// Generator of G2.
    pub prepared_h: PG::G2PreparedGadget,
    /// Generator of G1, times first monomial.
    pub prepared_beta_h: PG::G2PreparedGadget,
    /// Used for the shift powers associated with different degree bounds.
    pub prepared_degree_bounds_and_shift_powers: Option<Vec<(usize, Vec<PG::G1Gadget>)>>,
    np_phantom: PhantomData<NP>,
}

impl<CycleE, NP, PG> PreparedVerifierKeyGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    /// Find the appropriate shift for the degree bound.
    pub fn get_shift_power(&self, bound: usize) -> Option<Vec<PG::G1Gadget>> {
        self.prepared_degree_bounds_and_shift_powers
            .as_ref()
            .and_then(|v| {
                v.binary_search_by(|(d, _)| d.cmp(&bound))
                    .ok()
                    .map(|i| v[i].1.clone())
            })
    }
}

impl<CycleE, NP, PG>
    PrepareGadget<VerifierKeyGadget<CycleE, NP, PG>, <CycleE::E1 as PairingEngine>::Fr>
    for PreparedVerifierKeyGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn prepare<CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        unprepared: &VerifierKeyGadget<CycleE, NP, PG>,
    ) -> Result<Self, SynthesisError> {
        let supported_bits = <CycleE::E2 as PairingEngine>::Fr::size_in_bits();
        let mut prepared_g = Vec::<PG::G1Gadget>::new();

        let mut g: PG::G1Gadget = unprepared.g.clone();
        for i in 0..supported_bits {
            prepared_g.push(g.clone());
            g.double_in_place(&mut cs.ns(|| format!("g_double#{}", i)))?;
        }

        let prepared_h = PG::prepare_g2(&mut cs.ns(|| "h"), &unprepared.h)?;
        let prepared_beta_h = PG::prepare_g2(&mut cs.ns(|| "beta_h"), &unprepared.beta_h)?;

        let prepared_degree_bounds_and_shift_powers =
            if unprepared.degree_bounds_and_shift_powers.is_some() {
                let mut res = Vec::<(usize, Vec<PG::G1Gadget>)>::new();

                for (i, (d, shift_power)) in unprepared
                    .degree_bounds_and_shift_powers
                    .as_ref()
                    .unwrap()
                    .iter()
                    .enumerate()
                {
                    let mut prepared_shift_gadgets = Vec::<PG::G1Gadget>::new();

                    let mut cur: PG::G1Gadget = shift_power.clone();
                    for j in 0..supported_bits {
                        prepared_shift_gadgets.push(cur.clone());
                        cur.double_in_place(
                            &mut cs.ns(|| format!("shift_powers_double#{}, {}", i, j)),
                        )?;
                    }

                    res.push((d.clone(), prepared_shift_gadgets));
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
            np_phantom: PhantomData,
        })
    }
}

impl<CycleE, NP, PG> Clone for PreparedVerifierKeyGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn clone(&self) -> Self {
        PreparedVerifierKeyGadget {
            prepared_g: self.prepared_g.clone(),
            prepared_h: self.prepared_h.clone(),
            prepared_beta_h: self.prepared_beta_h.clone(),
            prepared_degree_bounds_and_shift_powers: self
                .prepared_degree_bounds_and_shift_powers
                .clone(),
            np_phantom: PhantomData,
        }
    }
}

impl<CycleE, NP, PG> AllocGadget<PreparedVerifierKey<CycleE::E2>, <CycleE::E1 as PairingEngine>::Fr>
    for PreparedVerifierKeyGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn alloc_constant<T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        t: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<PreparedVerifierKey<CycleE::E2>>,
    {
        let obj = t.borrow();

        let mut prepared_g = Vec::<PG::G1Gadget>::new();
        for (i, g) in obj.prepared_vk.prepared_g.iter().enumerate() {
            prepared_g.push(<PG::G1Gadget as AllocGadget<
                <CycleE::E2 as PairingEngine>::G1Projective,
                <CycleE::E1 as PairingEngine>::Fr,
            >>::alloc_constant(
                &mut cs.ns(|| format!("g#{}", i)),
                <CycleE::E2 as PairingEngine>::G1Projective::from(g.clone()),
            )?);
        }

        let prepared_h =
            PG::G2PreparedGadget::alloc_constant(&mut cs.ns(|| "h"), &obj.prepared_vk.prepared_h)?;
        let prepared_beta_h = PG::G2PreparedGadget::alloc_constant(
            &mut cs.ns(|| "beta_h"),
            &obj.prepared_vk.prepared_beta_h,
        )?;

        let prepared_degree_bounds_and_shift_powers = if obj
            .prepared_degree_bounds_and_shift_powers
            .is_some()
        {
            let mut res = Vec::<(usize, Vec<PG::G1Gadget>)>::new();

            for (i, (d, shift_power_elems)) in obj
                .prepared_degree_bounds_and_shift_powers
                .as_ref()
                .unwrap()
                .iter()
                .enumerate()
            {
                let mut gadgets = Vec::<PG::G1Gadget>::new();
                for (j, shift_power_elem) in shift_power_elems.iter().enumerate() {
                    gadgets.push(<PG::G1Gadget as AllocGadget<
                        <CycleE::E2 as PairingEngine>::G1Projective,
                        <CycleE::E1 as PairingEngine>::Fr,
                    >>::alloc_constant(
                        &mut cs.ns(|| format!("shift_power#{},{}", i, j)),
                        <CycleE::E2 as PairingEngine>::G1Projective::from(shift_power_elem.clone()),
                    )?);
                }

                res.push((d.clone(), gadgets));
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
            np_phantom: PhantomData,
        })
    }

    fn alloc<F, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        _cs: CS,
        _f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<PreparedVerifierKey<CycleE::E2>>,
    {
        unimplemented!()
    }

    fn alloc_input<F, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        _cs: CS,
        _f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<PreparedVerifierKey<CycleE::E2>>,
    {
        unimplemented!()
    }
}

/// Gadget for an optionally hiding Marlin-KZG10 commitment.
pub struct CommitmentGadget<
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
> where
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    comm: PG::G1Gadget,
    shifted_comm: Option<PG::G1Gadget>,
    np_phantom: PhantomData<NP>,
}

impl<CycleE, NP, PG> Clone for CommitmentGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn clone(&self) -> Self {
        CommitmentGadget {
            comm: self.comm.clone(),
            shifted_comm: self.shifted_comm.clone(),
            np_phantom: PhantomData,
        }
    }
}

impl<CycleE, NP, PG> AllocGadget<Commitment<CycleE::E2>, <CycleE::E1 as PairingEngine>::Fr>
    for CommitmentGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn alloc_constant<T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        _cs: CS,
        _val: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<Commitment<CycleE::E2>>,
    {
        unimplemented!()
    }

    #[inline]
    fn alloc<FN, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Commitment<CycleE::E2>>,
    {
        value_gen().and_then(|commitment| {
            let commitment = (*commitment.borrow()).clone();
            let comm = commitment.comm;
            let comm_gadget =
                PG::G1Gadget::alloc(cs.ns(|| "commitment"), || Ok(comm.0.into_projective()))?;

            let shifted_comm = commitment.shifted_comm;
            let shifted_comm_gadget = if let Some(shifted_comm) = shifted_comm {
                Some(PG::G1Gadget::alloc(cs.ns(|| "shifted commitment"), || {
                    Ok(shifted_comm.0.into_projective())
                })?)
            } else {
                None
            };

            Ok(Self {
                comm: comm_gadget,
                shifted_comm: shifted_comm_gadget,
                np_phantom: PhantomData,
            })
        })
    }

    #[inline]
    fn alloc_input<FN, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Commitment<CycleE::E2>>,
    {
        value_gen().and_then(|commitment| {
            let commitment = (*commitment.borrow()).clone();
            let comm = commitment.comm;
            let comm_gadget =
                PG::G1Gadget::alloc(cs.ns(|| "commitment"), || Ok(comm.0.into_projective()))?;

            let shifted_comm = commitment.shifted_comm;
            let shifted_comm_gadget = if let Some(shifted_comm) = shifted_comm {
                Some(PG::G1Gadget::alloc(cs.ns(|| "shifted commitment"), || {
                    Ok(shifted_comm.0.into_projective())
                })?)
            } else {
                None
            };

            Ok(Self {
                comm: comm_gadget,
                shifted_comm: shifted_comm_gadget,
                np_phantom: PhantomData,
            })
        })
    }
}

impl<CycleE, NP, PG> ToConstraintFieldGadget<<CycleE::E1 as PairingEngine>::Fr>
    for CommitmentGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn to_field_gadgets<CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<FpGadget<<CycleE::E1 as PairingEngine>::Fr>>, SynthesisError> {
        let mut res = Vec::new();

        let mut comm_gadget = self.comm.to_field_gadgets(&mut cs.ns(|| "comm"))?;

        res.append(&mut comm_gadget);

        if self.shifted_comm.as_ref().is_some() {
            let mut shifted_comm_gadget = self
                .shifted_comm
                .as_ref()
                .unwrap()
                .to_field_gadgets(&mut cs.ns(|| "shift_comm"))?;
            res.append(&mut shifted_comm_gadget);
        }

        Ok(res)
    }
}

impl<CycleE, NP, PG> ToBytesGadget<<CycleE::E1 as PairingEngine>::Fr>
    for CommitmentGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn to_bytes<CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let zero_shifted_comm = PG::G1Gadget::zero(&mut cs.ns(|| "zero_shifted_comm"))?;

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.comm.to_bytes(&mut cs.ns(|| "comm to bytes"))?);

        let shifted_comm = self.shifted_comm.clone().unwrap_or(zero_shifted_comm);
        bytes.extend_from_slice(&shifted_comm.to_bytes(&mut cs.ns(|| "shifted_comm to bytes"))?);
        Ok(bytes)
    }
}

/// Prepared gadget for an optionally hiding Marlin-KZG10 commitment.
/// shifted_comm is not prepared, due to the specific use case.
pub struct PreparedCommitmentGadget<
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
> where
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    prepared_comm: Vec<PG::G1Gadget>,
    shifted_comm: Option<PG::G1Gadget>,
    np_phantom: PhantomData<NP>,
}

impl<CycleE, NP, PG> Clone for PreparedCommitmentGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn clone(&self) -> Self {
        PreparedCommitmentGadget {
            prepared_comm: self.prepared_comm.clone(),
            shifted_comm: self.shifted_comm.clone(),
            np_phantom: PhantomData,
        }
    }
}

impl<CycleE, NP, PG>
    PrepareGadget<CommitmentGadget<CycleE, NP, PG>, <CycleE::E1 as PairingEngine>::Fr>
    for PreparedCommitmentGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn prepare<CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        unprepared: &CommitmentGadget<CycleE, NP, PG>,
    ) -> Result<Self, SynthesisError> {
        let mut prepared_comm = Vec::<PG::G1Gadget>::new();
        let supported_bits = <CycleE::E2 as PairingEngine>::Fr::size_in_bits();

        let mut cur: PG::G1Gadget = unprepared.comm.clone();
        for i in 0..supported_bits {
            prepared_comm.push(cur.clone());
            cur.double_in_place(&mut cs.ns(|| format!("double_in_place#{}", i)))?;
        }

        Ok(Self {
            prepared_comm,
            shifted_comm: unprepared.shifted_comm.clone(),
            np_phantom: PhantomData,
        })
    }
}

impl<CycleE, NP, PG> AllocGadget<PreparedCommitment<CycleE::E2>, <CycleE::E1 as PairingEngine>::Fr>
    for PreparedCommitmentGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn alloc_constant<T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        t: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<PreparedCommitment<CycleE::E2>>,
    {
        let obj = t.borrow();

        let mut prepared_comm = Vec::<PG::G1Gadget>::new();

        for (i, comm_elem) in obj.prepared_comm.0.iter().enumerate() {
            prepared_comm.push(<PG::G1Gadget as AllocGadget<
                <CycleE::E2 as PairingEngine>::G1Projective,
                <CycleE::E1 as PairingEngine>::Fr,
            >>::alloc_constant(
                &mut cs.ns(|| format!("comm_elem#{}", i)),
                <CycleE::E2 as PairingEngine>::G1Projective::from(comm_elem.clone()),
            )?);
        }

        let shifted_comm = if obj.shifted_comm.is_some() {
            Some(<PG::G1Gadget as AllocGadget<
                <CycleE::E2 as PairingEngine>::G1Projective,
                <CycleE::E1 as PairingEngine>::Fr,
            >>::alloc_constant(
                &mut cs.ns(|| "shifted_comm"),
                <CycleE::E2 as PairingEngine>::G1Projective::from(
                    obj.shifted_comm.unwrap().0.clone(),
                ),
            )?)
        } else {
            None
        };

        Ok(Self {
            prepared_comm,
            shifted_comm,
            np_phantom: PhantomData,
        })
    }

    fn alloc<F, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        _cs: CS,
        _f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<PreparedCommitment<CycleE::E2>>,
    {
        unimplemented!()
    }

    fn alloc_input<F, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        _cs: CS,
        _f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<PreparedCommitment<CycleE::E2>>,
    {
        unimplemented!()
    }
}

/// Gadget for a Marlin-KZG10 commitment, with a string label and degree bound.
pub struct LabeledCommitmentGadget<
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
> where
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    /// A text label for the commitment.
    pub label: String,
    /// The plain commitment.
    pub commitment: CommitmentGadget<CycleE, NP, PG>,
    /// Optionally, a bound on the polynomial degree.
    pub degree_bound: Option<UInt32>,
}

impl<CycleE, NP, PG> Clone for LabeledCommitmentGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn clone(&self) -> Self {
        LabeledCommitmentGadget {
            label: self.label.clone(),
            commitment: self.commitment.clone(),
            degree_bound: self.degree_bound.clone(),
        }
    }
}

impl<CycleE, NP, PG>
    AllocGadget<LabeledCommitment<Commitment<CycleE::E2>>, <CycleE::E1 as PairingEngine>::Fr>
    for LabeledCommitmentGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn alloc_constant<T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        _cs: CS,
        _val: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<LabeledCommitment<Commitment<CycleE::E2>>>,
    {
        unimplemented!()
    }

    #[inline]
    fn alloc<FN, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<LabeledCommitment<Commitment<CycleE::E2>>>,
    {
        value_gen().and_then(|labeled_commitment| {
            let labeled_commitment = labeled_commitment.borrow().clone();
            let label = labeled_commitment.label().to_string();
            let commitment = labeled_commitment.commitment();
            let degree_bound = labeled_commitment.degree_bound();

            let commitment =
                CommitmentGadget::alloc(cs.ns(|| format!("commitment {}", label)), || {
                    Ok(commitment)
                })?;

            let degree_bound = if degree_bound.is_some() {
                UInt32::alloc(
                    cs.ns(|| format!("degree_bound for {}", label)),
                    Some(degree_bound.unwrap() as u32),
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

    #[inline]
    fn alloc_input<FN, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<LabeledCommitment<Commitment<CycleE::E2>>>,
    {
        value_gen().and_then(|labeled_commitment| {
            let labeled_commitment = labeled_commitment.borrow().clone();
            let label = labeled_commitment.label().to_string();
            let commitment = labeled_commitment.commitment();
            let degree_bound = labeled_commitment.degree_bound();

            let commitment =
                CommitmentGadget::alloc(cs.ns(|| format!("commitment {}", label)), || {
                    Ok(commitment)
                })?;

            let degree_bound = if degree_bound.is_some() {
                UInt32::alloc(
                    cs.ns(|| format!("degree_bound for {}", label)),
                    Some(degree_bound.unwrap() as u32),
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

/// Gadget for a Marlin-KZG10 commitment, with a string label and degree bound.
pub struct PreparedLabeledCommitmentGadget<
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
> where
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    /// A text label for the commitment.
    pub label: String,
    /// The plain commitment.
    pub prepared_commitment: PreparedCommitmentGadget<CycleE, NP, PG>,
    /// Optionally, a bound on the polynomial degree.
    pub degree_bound: Option<UInt32>,
    np_phantom: PhantomData<NP>,
}

impl<CycleE, NP, PG> Clone for PreparedLabeledCommitmentGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn clone(&self) -> Self {
        PreparedLabeledCommitmentGadget {
            label: self.label.clone(),
            prepared_commitment: self.prepared_commitment.clone(),
            degree_bound: self.degree_bound.clone(),
            np_phantom: PhantomData,
        }
    }
}

impl<CycleE, NP, PG>
    PrepareGadget<LabeledCommitmentGadget<CycleE, NP, PG>, <CycleE::E1 as PairingEngine>::Fr>
    for PreparedLabeledCommitmentGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn prepare<CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        unprepared: &LabeledCommitmentGadget<CycleE, NP, PG>,
    ) -> Result<Self, SynthesisError> {
        let prepared_commitment = PreparedCommitmentGadget::prepare(
            cs.ns(|| "prepare unlabeled"),
            &unprepared.commitment,
        )?;

        Ok(Self {
            label: unprepared.label.clone(),
            prepared_commitment,
            degree_bound: unprepared.degree_bound.clone(),
            np_phantom: PhantomData,
        })
    }
}

/// Gadget for a Marlin-KZG10 proof.
pub struct ProofGadget<
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
> where
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    /// The commitment to the witness polynomial.
    pub w: PG::G1Gadget,
    /// The evaluation of the random hiding polynomial.
    pub random_v: Option<NonNativeFieldGadget<NP>>,
}

impl<CycleE, NP, PG> Clone for ProofGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn clone(&self) -> Self {
        ProofGadget {
            w: self.w.clone(),
            random_v: self.random_v.clone(),
        }
    }
}

impl<CycleE, NP, PG> AllocGadget<Proof<CycleE::E2>, <CycleE::E1 as PairingEngine>::Fr>
    for ProofGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn alloc_constant<T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        _cs: CS,
        _val: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<Proof<CycleE::E2>>,
    {
        unimplemented!()
    }

    #[inline]
    fn alloc<FN, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Proof<CycleE::E2>>,
    {
        value_gen().and_then(|proof| {
            let Proof { w, random_v } = proof.borrow().clone();
            let w = PG::G1Gadget::alloc(cs.ns(|| "w"), || Ok(w.into_projective()))?;

            let random_v = match random_v {
                None => None,
                Some(random_v_inner) => {
                    Some(NonNativeFieldGadget::alloc(cs.ns(|| "random_v"), || {
                        Ok(random_v_inner)
                    })?)
                }
            };

            Ok(Self { w, random_v })
        })
    }

    #[inline]
    fn alloc_input<FN, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Proof<CycleE::E2>>,
    {
        value_gen().and_then(|proof| {
            let Proof { w, random_v } = proof.borrow().clone();
            let w = PG::G1Gadget::alloc_input(cs.ns(|| "w"), || Ok(w.into_projective()))?;

            let random_v = match random_v {
                None => None,
                Some(random_v_inner) => {
                    Some(NonNativeFieldGadget::alloc(cs.ns(|| "random_v"), || {
                        Ok(random_v_inner)
                    })?)
                }
            };

            Ok(Self { w, random_v })
        })
    }
}

/// An allocated version of `LinearCombination`.
pub struct LinearCombinationGadget<
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
> where
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    /// The label.
    pub label: String,
    /// The linear combination of `(coeff, poly_label)` pairs.
    pub terms: Vec<(NonNativeFieldGadget<NP>, LCTerm)>,

    _cycle_engine: PhantomData<CycleE>,
    _pairing_gadget: PhantomData<PG>,
}

impl<CycleE, NP, PG> Clone for LinearCombinationGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn clone(&self) -> Self {
        LinearCombinationGadget {
            label: self.label.clone(),
            terms: self.terms.clone(),
            _cycle_engine: PhantomData,
            _pairing_gadget: PhantomData,
        }
    }
}

impl<CycleE, NP, PG>
    AllocGadget<
        LinearCombination<<CycleE::E2 as PairingEngine>::Fr>,
        <CycleE::E1 as PairingEngine>::Fr,
    > for LinearCombinationGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn alloc_constant<T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        val: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<LinearCombination<<CycleE::E2 as PairingEngine>::Fr>>,
    {
        let LinearCombination { label, terms } = val.borrow().clone();

        let new_terms: Vec<(NonNativeFieldGadget<NP>, LCTerm)> = terms
            .iter()
            .enumerate()
            .map(|(i, term)| {
                let (f, lc_term) = term;

                let fg = NonNativeFieldGadget::alloc_constant(cs.ns(|| format!("term {}", i)), f)
                    .unwrap();

                (fg, lc_term.clone())
            })
            .collect();

        Ok(Self {
            label: label.clone(),
            terms: new_terms,
            _cycle_engine: PhantomData,
            _pairing_gadget: PhantomData,
        })
    }

    #[inline]
    fn alloc<FN, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<LinearCombination<<CycleE::E2 as PairingEngine>::Fr>>,
    {
        value_gen().and_then(|lc| {
            let LinearCombination { label, terms } = lc.borrow().clone();

            let new_terms: Vec<(NonNativeFieldGadget<NP>, LCTerm)> = terms
                .iter()
                .enumerate()
                .map(|(i, term)| {
                    let (f, lc_term) = term;

                    let fg = NonNativeFieldGadget::alloc(cs.ns(|| format!("term {}", i)), || Ok(f))
                        .unwrap();

                    (fg, lc_term.clone())
                })
                .collect();

            Ok(Self {
                label: label.clone(),
                terms: new_terms,
                _cycle_engine: PhantomData,
                _pairing_gadget: PhantomData,
            })
        })
    }

    #[inline]
    fn alloc_input<FN, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<LinearCombination<<CycleE::E2 as PairingEngine>::Fr>>,
    {
        value_gen().and_then(|lc| {
            let LinearCombination { label, terms } = lc.borrow().clone();

            let new_terms: Vec<(NonNativeFieldGadget<NP>, LCTerm)> = terms
                .iter()
                .enumerate()
                .map(|(i, term)| {
                    let (f, lc_term) = term;

                    let fg = NonNativeFieldGadget::alloc(cs.ns(|| format!("term {}", i)), || Ok(f))
                        .unwrap();

                    (fg, lc_term.clone())
                })
                .collect();

            Ok(Self {
                label: label.clone(),
                terms: new_terms,
                _cycle_engine: PhantomData,
                _pairing_gadget: PhantomData,
            })
        })
    }
}

/// An allocated version of `BatchLCProof`.
pub struct BatchLCProofGadget<
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
> where
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    /// Evaluation proofs.
    pub proofs: Vec<ProofGadget<CycleE, NP, PG>>,
    /// Evaluations required to verify the proof.
    pub evals: Option<Vec<NonNativeFieldGadget<NP>>>,
}

impl<CycleE, NP, PG> Clone for BatchLCProofGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn clone(&self) -> Self {
        BatchLCProofGadget {
            proofs: self.proofs.clone(),
            evals: self.evals.clone(),
        }
    }
}

impl<CycleE, NP, PG>
    AllocGadget<
        BatchLCProof<<CycleE::E2 as PairingEngine>::Fr, MarlinKZG10<CycleE::E2>>,
        <CycleE::E1 as PairingEngine>::Fr,
    > for BatchLCProofGadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    fn alloc_constant<T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        _cs: CS,
        _val: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<BatchLCProof<<CycleE::E2 as PairingEngine>::Fr, MarlinKZG10<CycleE::E2>>>,
    {
        unimplemented!()
    }

    #[inline]
    fn alloc<FN, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<BatchLCProof<<CycleE::E2 as PairingEngine>::Fr, MarlinKZG10<CycleE::E2>>>,
    {
        value_gen().and_then(|proof| {
            let BatchLCProof { proof, evals } = proof.borrow().clone();

            let proofs: Vec<Proof<_>> = proof.clone().into();
            let proofs: Vec<ProofGadget<CycleE, NP, PG>> = proofs
                .iter()
                .enumerate()
                .map(|(i, p)| {
                    ProofGadget::alloc(cs.ns(|| format!("proof {}", i)), || Ok(p)).unwrap()
                })
                .collect();

            let evals: Option<Vec<NonNativeFieldGadget<NP>>> = match evals {
                None => None,
                Some(evals_inner) => Some(
                    evals_inner
                        .iter()
                        .enumerate()
                        .map(|(i, e)| {
                            NonNativeFieldGadget::alloc(
                                cs.ns(|| format!("evaluation {}", i)),
                                || Ok(e),
                            )
                            .unwrap()
                        })
                        .collect(),
                ),
            };

            Ok(Self { proofs, evals })
        })
    }

    #[inline]
    fn alloc_input<FN, T, CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<BatchLCProof<<CycleE::E2 as PairingEngine>::Fr, MarlinKZG10<CycleE::E2>>>,
    {
        value_gen().and_then(|proof| {
            let BatchLCProof { proof, evals } = proof.borrow().clone();

            let proofs: Vec<Proof<_>> = proof.clone().into();
            let proofs: Vec<ProofGadget<CycleE, NP, PG>> = proofs
                .iter()
                .enumerate()
                .map(|(i, p)| {
                    ProofGadget::alloc(cs.ns(|| format!("proof {}", i)), || Ok(p)).unwrap()
                })
                .collect();

            let evals: Option<Vec<NonNativeFieldGadget<NP>>> = match evals {
                None => None,
                Some(evals_inner) => Some(
                    evals_inner
                        .iter()
                        .enumerate()
                        .map(|(i, e)| {
                            NonNativeFieldGadget::alloc(
                                cs.ns(|| format!("evaluation {}", i)),
                                || Ok(e),
                            )
                            .unwrap()
                        })
                        .collect(),
                ),
            };

            Ok(Self { proofs, evals })
        })
    }
}

/// Gadget for the Marlin-KZG10 polynomial commitment verifier.
pub struct MarlinKZG10Gadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    _cycle_engine: PhantomData<CycleE>,
    _nonnative_params_gadget: PhantomData<NP>,
    _pairing_gadget: PhantomData<PG>,
}

impl<CycleE, NP, PG>
    PCCheckGadget<
        <CycleE::E2 as PairingEngine>::Fr,
        MarlinKZG10<CycleE::E2>,
        NP,
        <CycleE::E1 as PairingEngine>::Fr,
    > for MarlinKZG10Gadget<CycleE, NP, PG>
where
    CycleE: CycleEngine,
    NP: NonNativeFieldParams<
        BaseField = <CycleE::E1 as PairingEngine>::Fr,
        TargetField = <CycleE::E2 as PairingEngine>::Fr,
    >,
    PG: PairingGadget<CycleE::E2, <CycleE::E1 as PairingEngine>::Fr>,
    <CycleE::E2 as PairingEngine>::G1Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
    <CycleE::E2 as PairingEngine>::G2Projective: MulAssign<<CycleE::E1 as PairingEngine>::Fq>,
{
    type VerifierKeyGadget = VerifierKeyGadget<CycleE, NP, PG>;
    type PreparedVerifierKeyGadget = PreparedVerifierKeyGadget<CycleE, NP, PG>;
    type CommitmentGadget = CommitmentGadget<CycleE, NP, PG>;
    type PreparedCommitmentGadget = PreparedCommitmentGadget<CycleE, NP, PG>;
    type LabeledCommitmentGadget = LabeledCommitmentGadget<CycleE, NP, PG>;
    type PreparedLabeledCommitmentGadget = PreparedLabeledCommitmentGadget<CycleE, NP, PG>;
    type ProofGadget = ProofGadget<CycleE, NP, PG>;
    type LinearCombinationGadget = LinearCombinationGadget<CycleE, NP, PG>;
    type BatchLCProofGadget = BatchLCProofGadget<CycleE, NP, PG>;

    fn batch_check_evaluations<CS>(
        mut cs: CS,
        verification_key: &Self::VerifierKeyGadget,
        commitments: &Vec<Self::LabeledCommitmentGadget>,
        query_set: &QuerySetGadget<NP>,
        evaluations: &EvaluationsGadget<NP>,
        proofs: &Vec<Self::ProofGadget>,
        opening_challenges: &Vec<NonNativeFieldGadget<NP>>,
        opening_challenges_bits: &Vec<Vec<Boolean>>,
        batching_rands: &Vec<NonNativeFieldGadget<NP>>,
        batching_rands_bits: &Vec<Vec<Boolean>>,
    ) -> Result<(), SynthesisError>
    where
        CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>,
    {
        let initializer =
            <<CycleE as CycleEngine>::E2 as PairingEngine>::G1Projective::prime_subgroup_generator(
            );
        let initializer_gadget = PG::G1Gadget::alloc(
            cs.ns(|| "batch_check_evaluations: the initializer used to combat incomplete group operations"),
            || Ok(initializer),
        )?;

        let mut batching_rands = batching_rands.clone();
        let mut batching_rands_bits = batching_rands_bits.clone();

        let commitments: BTreeMap<_, _> =
            commitments.iter().map(|c| (c.label.clone(), c)).collect();
        let mut query_to_labels_map = BTreeMap::new();

        for (label, point) in query_set.iter() {
            let labels = query_to_labels_map.entry(point).or_insert(BTreeSet::new());
            labels.insert(label);
        }

        let mut combined_comms = Vec::new();
        let mut combined_queries = Vec::new();
        let mut combined_evals = Vec::new();
        for (c, (query, labels)) in query_to_labels_map.into_iter().enumerate() {
            let mut comms_to_combine: Vec<Self::LabeledCommitmentGadget> = Vec::new();
            let mut values_to_combine = Vec::new();
            for label in labels.into_iter() {
                let commitment = commitments.get(label).unwrap().clone();
                let degree_bound = commitment.degree_bound.clone();
                assert_eq!(
                    degree_bound.is_some(),
                    commitment.commitment.shifted_comm.is_some()
                );

                let v_i = evaluations.get(&(label.clone(), query.clone())).unwrap();

                comms_to_combine.push(commitment.clone());
                values_to_combine.push(v_i.clone());
            }

            let mut cs = cs.ns(|| format!("process batch {}", c));

            // begin accumulate_commitments_and_values
            let mut combined_comm = initializer_gadget.clone();
            let mut combined_value =
                NonNativeFieldMulResultGadget::<NP>::zero(cs.ns(|| "zero for combined_value"))?;

            let mut opening_challenges_counter = 0;

            for (i, (labeled_commitment, value)) in comms_to_combine
                .into_iter()
                .zip(values_to_combine.iter())
                .enumerate()
            {
                let challenge = opening_challenges[opening_challenges_counter].clone();
                let challenge_bits = opening_challenges_bits[opening_challenges_counter].clone();
                opening_challenges_counter += 1;

                let mut cs = cs.ns(|| format!("process commitment {}", i));

                let LabeledCommitmentGadget {
                    commitment,
                    degree_bound,
                    ..
                } = labeled_commitment;
                let CommitmentGadget { shifted_comm, .. } = commitment;

                combined_comm = commitment.comm.mul_bits(
                    cs.ns(|| "comm * challenge"),
                    &combined_comm,
                    challenge_bits.iter().rev(),
                )?;

                let value_times_challenge_unreduced =
                    value.mul_without_reduce(&mut cs.ns(|| "value * challenge"), &challenge)?;
                combined_value = combined_value.add(
                    cs.ns(|| "combined_value += value * challenge"),
                    &value_times_challenge_unreduced,
                )?;

                if let Some(degree_bound) = degree_bound {
                    let challenge_shifted_bits =
                        opening_challenges_bits[opening_challenges_counter].clone();
                    opening_challenges_counter += 1;

                    let shifted_comm = shifted_comm.as_ref().unwrap();

                    let value_bits = value.to_bits(cs.ns(|| "value to bits"))?;
                    let shift_power = verification_key
                        .get_shift_power(degree_bound.value.unwrap() as usize)
                        .unwrap();

                    let mut adjusted_comm = shift_power.mul_bits(
                        cs.ns(|| "adjusted_comm"),
                        &initializer_gadget,
                        value_bits.iter().rev(),
                    )?;
                    adjusted_comm = adjusted_comm.sub(
                        cs.ns(|| "subtract rand for adjusted_comm"),
                        &initializer_gadget,
                    )?;

                    adjusted_comm = shifted_comm.sub(
                        cs.ns(|| "shifted_comm - shift_power * value"),
                        &adjusted_comm,
                    )?;

                    let mut final_adjusted_comm = adjusted_comm.mul_bits(
                        cs.ns(|| "final_adjusted_comm"),
                        &initializer_gadget,
                        challenge_shifted_bits.iter().rev(),
                    )?;
                    final_adjusted_comm = final_adjusted_comm.sub(
                        cs.ns(|| "subtract rand for final_adjusted_comm"),
                        &initializer_gadget,
                    )?;

                    combined_comm = combined_comm.add(
                        cs.ns(|| "combined_comm += final_adjusted_comm"),
                        &final_adjusted_comm,
                    )?;
                }
            }

            combined_comm =
                combined_comm.sub(cs.ns(|| "subtract rand out"), &initializer_gadget)?;

            combined_comms.push(combined_comm);
            combined_queries.push(query);
            combined_evals.push(combined_value);
        }

        {
            let mut cs = cs.ns(|| "batch check");

            let mut total_c = initializer_gadget.clone();
            let mut total_w = initializer_gadget.clone();

            let mut g_multiplier =
                NonNativeFieldMulResultGadget::<NP>::zero(cs.ns(|| "zero for g_multiplier"))?;
            for (i, (((c, z), v), proof)) in combined_comms
                .iter()
                .zip(combined_queries)
                .zip(combined_evals)
                .zip(proofs)
                .enumerate()
            {
                let mut cs = cs.ns(|| format!("process combined commitment {}", i));

                let z_bits = z.to_bits(cs.ns(|| "z to bits"))?;

                let mut w_times_z = proof.w.mul_bits(
                    cs.ns(|| "w_times_z"),
                    &initializer_gadget,
                    z_bits.iter().rev(),
                )?;
                w_times_z =
                    w_times_z.sub(cs.ns(|| "subtract rand for w_times_z"), &initializer_gadget)?;

                let c = c.add(cs.ns(|| "c += w * z"), &w_times_z)?;

                let randomizer = batching_rands.remove(0);
                let randomizer_bits = batching_rands_bits.remove(0);

                let v_reduced = v.reduce(&mut cs.ns(|| "reduced v"))?;
                let randomizer_times_v =
                    randomizer.mul_without_reduce(&mut cs.ns(|| "randomizer * v"), &v_reduced)?;

                g_multiplier = g_multiplier.add(
                    cs.ns(|| "g_mult += randomizer times v"),
                    &randomizer_times_v,
                )?;

                let mut c_times_randomizer = c.mul_bits(
                    cs.ns(|| "c times randomizer"),
                    &initializer_gadget,
                    randomizer_bits.iter().rev(),
                )?;
                c_times_randomizer = c_times_randomizer.sub(
                    cs.ns(|| "subtract rand for c times randomizer"),
                    &initializer_gadget,
                )?;

                let mut w_times_randomizer = proof.w.mul_bits(
                    cs.ns(|| "w times randomizer"),
                    &initializer_gadget,
                    randomizer_bits.iter().rev(),
                )?;
                w_times_randomizer = w_times_randomizer.sub(
                    cs.ns(|| "subtract rand for w times randomizer"),
                    &initializer_gadget,
                )?;

                total_c = total_c.add(
                    cs.ns(|| format!("add to c for round {}", i)),
                    &c_times_randomizer,
                )?;

                total_w = total_w.add(
                    cs.ns(|| format!("add to w for round {}", i)),
                    &w_times_randomizer,
                )?;
            }

            let (prepared_total_w, prepared_beta_h, prepared_total_c, prepared_h) = {
                let mut cs = cs.ns(|| "prepare for pairing check");

                let g_multiplier_reduced = g_multiplier.reduce(cs.ns(|| "reduce g_multiplier"))?;
                let g_multiplier_bits =
                    g_multiplier_reduced.to_bits(cs.ns(|| "g_multiplier to bits"))?;

                let mut g_times_mul = verification_key.g.mul_bits(
                    cs.ns(|| "g * g_multiplier"),
                    &initializer_gadget,
                    g_multiplier_bits.iter().rev(),
                )?;
                g_times_mul = g_times_mul.sub(
                    cs.ns(|| "subtract rand for g * g_multiplier"),
                    &initializer_gadget,
                )?;

                total_c = total_c.sub(
                    cs.ns(|| "subtract off initial rand total_c"),
                    &initializer_gadget,
                )?;
                total_w = total_w.sub(
                    cs.ns(|| "subtract off initial rand total_w"),
                    &initializer_gadget,
                )?;

                total_c = total_c.sub(
                    cs.ns(|| "subtract g * g_multiplier from total_c"),
                    &g_times_mul,
                )?;

                total_w = total_w.negate(cs.ns(|| "negate total_w"))?;

                let prepared_total_w =
                    PG::prepare_g1(cs.ns(|| "prepare total_w for lhs"), &total_w)?;
                let prepared_beta_h =
                    PG::prepare_g2(cs.ns(|| "prepare beta_h"), &verification_key.beta_h)?;
                let prepared_total_c =
                    PG::prepare_g1(cs.ns(|| "prepare total_c for lhs"), &total_c)?;
                let prepared_h = PG::prepare_g2(cs.ns(|| "prepare h"), &verification_key.h)?;

                (
                    prepared_total_w,
                    prepared_beta_h,
                    prepared_total_c,
                    prepared_h,
                )
            };

            let lhs = PG::product_of_pairings(
                cs.ns(|| "lhs = pairing(total_w, beta_h) * pairing(total_c, h)"),
                &[prepared_total_w, prepared_total_c],
                &[prepared_beta_h, prepared_h],
            )?;

            let rhs = &PG::GTGadget::one(cs.ns(|| "one for RHS"))?;

            lhs.enforce_equal(cs.ns(|| "lhs = rhs"), rhs)?;
        }

        Ok(())
    }

    fn prepared_batch_check_evaluations<CS>(
        mut cs: CS,
        prepared_verification_key: &Self::PreparedVerifierKeyGadget,
        prepared_commitments: &Vec<Self::PreparedLabeledCommitmentGadget>,
        query_set: &QuerySetGadget<NP>,
        evaluations: &EvaluationsGadget<NP>,
        proofs: &Vec<Self::ProofGadget>,
        opening_challenges: &Vec<NonNativeFieldGadget<NP>>,
        opening_challenges_bits: &Vec<Vec<Boolean>>,
        batching_rands: &Vec<NonNativeFieldGadget<NP>>,
        batching_rands_bits: &Vec<Vec<Boolean>>,
    ) -> Result<(), SynthesisError>
    where
        CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>,
    {
        let initializer =
            <<CycleE as CycleEngine>::E2 as PairingEngine>::G1Projective::prime_subgroup_generator(
            );
        let initializer_gadget = PG::G1Gadget::alloc(
            cs.ns(|| "prepared_batch_check_evaluations: the initializer used to combat incomplete group operations"),
            || Ok(initializer),
        )?;

        let mut batching_rands = batching_rands.clone();
        let mut batching_rands_bits = batching_rands_bits.clone();

        let commitments: BTreeMap<_, _> = prepared_commitments
            .iter()
            .map(|c| (c.label.clone(), c))
            .collect();
        let mut query_to_labels_map = BTreeMap::new();

        for (label, point) in query_set.iter() {
            let labels = query_to_labels_map.entry(point).or_insert(BTreeSet::new());
            labels.insert(label);
        }

        let mut combined_comms = Vec::new();
        let mut combined_queries = Vec::new();
        let mut combined_evals = Vec::new();
        for (c, (query, labels)) in query_to_labels_map.into_iter().enumerate() {
            let mut comms_to_combine: Vec<Self::PreparedLabeledCommitmentGadget> = Vec::new();
            let mut values_to_combine = Vec::new();
            for label in labels.into_iter() {
                let commitment = commitments.get(label).unwrap().clone();
                let degree_bound = commitment.degree_bound.clone();
                assert_eq!(
                    degree_bound.is_some(),
                    commitment.prepared_commitment.shifted_comm.is_some()
                );

                let v_i = evaluations.get(&(label.clone(), query.clone())).unwrap();

                comms_to_combine.push(commitment.clone());
                values_to_combine.push(v_i.clone());
            }

            let mut cs = cs.ns(|| format!("process batch {}", c));

            // begin accumulate_commitments_and_values
            let mut combined_comm = PG::G1Gadget::zero(cs.ns(|| "initialize combined commitment"))?;
            let mut combined_comm_is_zero = true;
            let mut combined_value =
                NonNativeFieldMulResultGadget::<NP>::zero(cs.ns(|| "zero for combined_value"))?;

            let mut opening_challenges_counter = 0;

            for (i, (labeled_commitment, value)) in comms_to_combine
                .into_iter()
                .zip(values_to_combine)
                .enumerate()
            {
                let challenge = opening_challenges[opening_challenges_counter].clone();
                let challenge_bits = opening_challenges_bits[opening_challenges_counter].clone();
                opening_challenges_counter += 1;

                let mut cs = cs.ns(|| format!("process commitment {}", i));

                let PreparedLabeledCommitmentGadget {
                    prepared_commitment,
                    degree_bound,
                    ..
                } = labeled_commitment;

                let PreparedCommitmentGadget { shifted_comm, .. } = prepared_commitment;

                let mut comm_times_challenge: PG::G1Gadget = initializer_gadget.clone();
                {
                    let cs = &mut cs.ns(|| "comm * challenge");
                    for (j, (bit, base_power)) in challenge_bits
                        .iter()
                        .rev()
                        .zip(&prepared_commitment.prepared_comm)
                        .enumerate()
                    {
                        let new_encoded = comm_times_challenge.add(
                            &mut cs.ns(|| format!("Add {}-th base power", j)),
                            &base_power,
                        )?;
                        comm_times_challenge = PG::G1Gadget::conditionally_select(
                            &mut cs.ns(|| format!("Conditional Select {}", j)),
                            bit,
                            &new_encoded,
                            &comm_times_challenge,
                        )?;
                    }
                }
                comm_times_challenge =
                    comm_times_challenge.sub(cs.ns(|| "subtract rand out"), &initializer_gadget)?;

                if combined_comm_is_zero {
                    combined_comm = comm_times_challenge;
                    combined_comm_is_zero = false;
                } else {
                    combined_comm = combined_comm.add(
                        cs.ns(|| "combined_comm += comm * challenge"),
                        &comm_times_challenge,
                    )?;
                }

                let value_times_challenge_unreduced =
                    value.mul_without_reduce(&mut cs.ns(|| "value * challenge"), &challenge)?;
                combined_value = combined_value.add(
                    cs.ns(|| "combined_value += value * challenge"),
                    &value_times_challenge_unreduced,
                )?;

                if let Some(degree_bound) = degree_bound {
                    let challenge_shifted_bits =
                        opening_challenges_bits[opening_challenges_counter].clone();
                    opening_challenges_counter += 1;

                    let shifted_comm = shifted_comm.as_ref().unwrap();

                    let value_bits = value.to_bits(cs.ns(|| "value_to_bits"))?;
                    let shift_power = prepared_verification_key
                        .get_shift_power(degree_bound.value.unwrap() as usize)
                        .unwrap();
                    //.ok_or(<MarlinKZG10<CycleE::E2> as PolynomialCommitment<<CycleE::E2 as PairingEngine>::Fr>>::Error::UnsupportedDegreeBound(0))?;
                    // maybetodo: use usize version of degree_bound sometime

                    let mut adjusted_comm = initializer_gadget.clone();
                    {
                        let cs = &mut cs.ns(|| "adjusted_comm");
                        for (j, (bit, base_power)) in
                            value_bits.iter().rev().zip(&shift_power).enumerate()
                        {
                            let new_encoded = adjusted_comm.add(
                                &mut cs.ns(|| format!("Add {}-th base power", j)),
                                &base_power,
                            )?;
                            adjusted_comm = PG::G1Gadget::conditionally_select(
                                &mut cs.ns(|| format!("Conditional Select {}", j)),
                                bit,
                                &new_encoded,
                                &adjusted_comm,
                            )?;
                        }
                    }
                    adjusted_comm = adjusted_comm.sub(
                        cs.ns(|| "subtract rand for adjusted_comm"),
                        &initializer_gadget,
                    )?;

                    adjusted_comm = shifted_comm.sub(
                        cs.ns(|| "shifted_comm - shift_power * value"),
                        &adjusted_comm,
                    )?;

                    let mut final_adjusted_comm = adjusted_comm.mul_bits(
                        cs.ns(|| "final_adjusted_comm"),
                        &initializer_gadget,
                        challenge_shifted_bits.iter().rev(),
                    )?;
                    final_adjusted_comm = final_adjusted_comm.sub(
                        cs.ns(|| "subtract rand for final_adjusted_comm"),
                        &initializer_gadget,
                    )?;

                    combined_comm = combined_comm.add(
                        cs.ns(|| "combined_comm += final_adjusted_comm"),
                        &final_adjusted_comm,
                    )?;
                }
            }

            combined_comms.push(combined_comm);
            combined_queries.push(query);
            combined_evals.push(combined_value);
        }

        {
            let mut cs = cs.ns(|| "batch check");

            let mut total_c = initializer_gadget.clone();
            let mut total_w = initializer_gadget.clone();

            let mut g_multiplier =
                NonNativeFieldMulResultGadget::<NP>::zero(cs.ns(|| "zero for g_multiplier"))?;
            for (i, (((c, z), v), proof)) in combined_comms
                .iter()
                .zip(combined_queries)
                .zip(combined_evals)
                .zip(proofs)
                .enumerate()
            {
                let mut cs = cs.ns(|| format!("process combined commitment {}", i));

                let z_bits = z.to_bits(cs.ns(|| "z to bits"))?;

                let mut w_times_z = proof.w.mul_bits(
                    cs.ns(|| "w_times_z"),
                    &initializer_gadget,
                    z_bits.iter().rev(),
                )?;
                w_times_z =
                    w_times_z.sub(cs.ns(|| "subtract rand for w_times_z"), &initializer_gadget)?;

                let c = c.add(cs.ns(|| "c += w * z"), &w_times_z)?;

                let randomizer = batching_rands.remove(0);
                let randomizer_bits = batching_rands_bits.remove(0);

                let v_reduced = v.reduce(&mut cs.ns(|| "reduced v"))?;
                let randomizer_times_v =
                    randomizer.mul_without_reduce(&mut cs.ns(|| "randomizer * v"), &v_reduced)?;

                g_multiplier = g_multiplier.add(
                    cs.ns(|| "g_mult += randomizer times v"),
                    &randomizer_times_v,
                )?;

                let mut c_times_randomizer = c.mul_bits(
                    cs.ns(|| "c times randomizer"),
                    &initializer_gadget,
                    randomizer_bits.iter().rev(),
                )?;
                c_times_randomizer = c_times_randomizer.sub(
                    cs.ns(|| "subtract rand for c times randomizer"),
                    &initializer_gadget,
                )?;

                let mut w_times_randomizer = proof.w.mul_bits(
                    cs.ns(|| "w times randomizer"),
                    &initializer_gadget,
                    randomizer_bits.iter().rev(),
                )?;
                w_times_randomizer = w_times_randomizer.sub(
                    cs.ns(|| "subtract rand for w times randomizer"),
                    &initializer_gadget,
                )?;

                total_c = total_c.add(
                    cs.ns(|| format!("add to c for round {}", i)),
                    &c_times_randomizer,
                )?;

                total_w = total_w.add(
                    cs.ns(|| format!("add to w for round {}", i)),
                    &w_times_randomizer,
                )?;
            }

            let (prepared_total_w, prepared_beta_h, prepared_total_c, prepared_h) = {
                let mut cs = cs.ns(|| "prepare for pairing check");

                let g_multiplier_reduced = g_multiplier.reduce(cs.ns(|| "reduce g_multiplier"))?;
                let g_multiplier_bits =
                    g_multiplier_reduced.to_bits(cs.ns(|| "g_multiplier to bits"))?;

                let mut g_times_mul = initializer_gadget.clone();
                {
                    let cs = &mut cs.ns(|| "g * g_multiplier");
                    for (j, (bit, base_power)) in g_multiplier_bits
                        .iter()
                        .rev()
                        .zip(&prepared_verification_key.prepared_g)
                        .enumerate()
                    {
                        let new_encoded = g_times_mul.add(
                            &mut cs.ns(|| format!("Add {}-th base power", j)),
                            &base_power,
                        )?;
                        g_times_mul = PG::G1Gadget::conditionally_select(
                            &mut cs.ns(|| format!("Conditional Select {}", j)),
                            bit,
                            &new_encoded,
                            &g_times_mul,
                        )?;
                    }
                }
                g_times_mul = g_times_mul.sub(
                    cs.ns(|| "subtract rand for g * g_multiplier"),
                    &initializer_gadget,
                )?;

                total_c = total_c.sub(
                    cs.ns(|| "subtract off initial rand total_c"),
                    &initializer_gadget,
                )?;
                total_w = total_w.sub(
                    cs.ns(|| "subtract off initial rand total_w"),
                    &initializer_gadget,
                )?;

                total_c = total_c.sub(
                    cs.ns(|| "subtract g * g_multiplier from total_c"),
                    &g_times_mul,
                )?;

                total_w = total_w.negate(cs.ns(|| "negate total_w"))?;

                let prepared_total_w =
                    PG::prepare_g1(cs.ns(|| "prepare total_w for lhs"), &total_w)?;
                let prepared_beta_h = prepared_verification_key.prepared_beta_h.clone();
                let prepared_total_c =
                    PG::prepare_g1(cs.ns(|| "prepare total_c for lhs"), &total_c)?;
                let prepared_h = prepared_verification_key.prepared_h.clone();

                (
                    prepared_total_w,
                    prepared_beta_h,
                    prepared_total_c,
                    prepared_h,
                )
            };

            let lhs = PG::product_of_pairings(
                cs.ns(|| "lhs = pairing(total_w, beta_h) * pairing(total_c, h)"),
                &[prepared_total_w, prepared_total_c],
                &[prepared_beta_h, prepared_h],
            )?;

            let rhs = &PG::GTGadget::one(cs.ns(|| "one for RHS"))?;

            lhs.enforce_equal(cs.ns(|| "lhs = rhs"), rhs)?;
        }

        Ok(())
    }

    fn prepared_check_combinations<CS>(
        mut cs: CS,
        prepared_verification_key: &Self::PreparedVerifierKeyGadget,
        linear_combinations: &Vec<Self::LinearCombinationGadget>,
        prepared_commitments: &Vec<Self::PreparedLabeledCommitmentGadget>,
        query_set: &QuerySetGadget<NP>,
        evaluations: &EvaluationsGadget<NP>,
        proof: &Self::BatchLCProofGadget,
        opening_challenges: &Vec<NonNativeFieldGadget<NP>>,
        opening_challenges_bits: &Vec<Vec<Boolean>>,
        batching_rands: &Vec<NonNativeFieldGadget<NP>>,
        batching_rands_bits: &Vec<Vec<Boolean>>,
    ) -> Result<(), SynthesisError>
    where
        CS: ConstraintSystem<<CycleE::E1 as PairingEngine>::Fr>,
    {
        let BatchLCProofGadget { proofs, .. } = proof;

        let initializer =
            <<CycleE as CycleEngine>::E2 as PairingEngine>::G1Projective::prime_subgroup_generator(
            );
        let initializer_gadget = PG::G1Gadget::alloc(
            cs.ns(|| "prepared_check_combinations: the initializer used to combat incomplete group operations"),
            || Ok(initializer),
        )?;

        let label_comm_map = prepared_commitments
            .iter()
            .map(|c| (c.label.clone(), c))
            .collect::<BTreeMap<_, _>>();

        let mut lc_commitments = Vec::new();
        let mut lc_info = Vec::new();
        let mut evaluations = evaluations.clone();

        for (lc_num, lc) in linear_combinations.iter().enumerate() {
            let lc_label = lc.label.clone();
            let num_polys = lc.terms.len();

            let mut degree_bound = None;
            let mut coeffs_and_comms = Vec::new();

            for (term_num, (coeff, label)) in lc.terms.iter().enumerate() {
                if label.is_one() {
                    for (&(ref label, _), ref mut eval) in evaluations.iter_mut() {
                        if label == &lc_label {
                            **eval = eval.sub(
                                cs.ns(|| format!("subtract eval for label {}", label)),
                                coeff,
                            )?;
                        }
                    }
                } else {
                    let label: &String = label.try_into().unwrap();
                    let &cur_comm = label_comm_map.get(label).unwrap();

                    if num_polys == 1 && cur_comm.degree_bound.is_some() {
                        let one = NonNativeFieldGadget::one(
                            cs.ns(|| format!("coeff = one for {}, {}", lc_num, term_num)),
                        )?;
                        coeff.enforce_equal(
                            cs.ns(|| {
                                format!(
                                    "coefficient must be one for degree-bounded equations {}, {}",
                                    lc_num, term_num
                                )
                            }),
                            &one,
                        )?;
                        degree_bound = cur_comm.degree_bound.clone();
                    }
                    coeffs_and_comms.push((coeff.clone(), cur_comm.prepared_commitment.clone()));
                }
            }

            let mut combined_comm: PG::G1Gadget = initializer_gadget.clone();
            let mut combined_shifted_comm: Option<PG::G1Gadget> = None;

            for (c, (coeff, comm)) in coeffs_and_comms.iter().enumerate() {
                let mut cs = cs.ns(|| format!("coeffs and comms, {}, round {}", lc_num, c));

                let coeff_bits = coeff.to_bits(cs.ns(|| "coeff to bits")).unwrap();

                let mut comm_times_coeff: PG::G1Gadget = initializer_gadget.clone();
                {
                    let cs = &mut cs.ns(|| "comm * challenge");
                    for (j, (bit, base_power)) in
                        coeff_bits.iter().rev().zip(&comm.prepared_comm).enumerate()
                    {
                        let new_encoded = comm_times_coeff.add(
                            &mut cs.ns(|| format!("Add {}-th base power", j)),
                            &base_power,
                        )?;
                        comm_times_coeff = PG::G1Gadget::conditionally_select(
                            &mut cs.ns(|| format!("Conditional Select {}", j)),
                            bit,
                            &new_encoded,
                            &comm_times_coeff,
                        )?;
                    }
                }
                comm_times_coeff = comm_times_coeff.sub(
                    cs.ns(|| format!("subtract rand out of comm * coeff for round {}", c)),
                    &initializer_gadget,
                )?;
                combined_comm = combined_comm.add(
                    cs.ns(|| format!("add in comm * coeff for round {}", c)),
                    &comm_times_coeff,
                )?;

                if let Some(shifted_comm) = &comm.shifted_comm {
                    let mut cur = shifted_comm.mul_bits(
                        cs.ns(|| "shifted_comm * coeff"),
                        &initializer_gadget,
                        coeff_bits.iter().rev(),
                    )?;
                    cur = cur.sub(
                        cs.ns(|| "subtract rand for shifted_comm * coeff"),
                        &initializer_gadget,
                    )?;

                    if let Some(combined_shifted_comm_inner) = combined_shifted_comm {
                        combined_shifted_comm = Some(
                            combined_shifted_comm_inner
                                .add(cs.ns(|| "combined_shifted_comm += cur"), &cur)?,
                        );
                    } else {
                        combined_shifted_comm = Some(cur);
                    }
                }
            }
            combined_comm = combined_comm.sub(
                cs.ns(|| format!("subtract rand out of combined_comm {}", lc_num)),
                &initializer_gadget,
            )?;
            lc_commitments.push((combined_comm, combined_shifted_comm));

            lc_info.push((lc_label, degree_bound));
        }

        //let comms = Self::normalize_commitments(lc_commitments);
        let comms = lc_commitments;

        let lc_commitments = lc_info
            .iter()
            .zip(comms)
            .enumerate()
            .map(|(i, ((label, d), (comm, shifted_comm)))| {
                let commitment = CommitmentGadget {
                    comm,
                    shifted_comm,
                    np_phantom: PhantomData,
                };
                let labeled_commitment = LabeledCommitmentGadget {
                    label: label.to_string(),
                    commitment,
                    degree_bound: d.clone(),
                };
                PreparedLabeledCommitmentGadget::prepare(
                    &mut cs.ns(|| format!("prepare commitment {}", i)),
                    &labeled_commitment,
                )
                .unwrap()
            })
            .collect::<Vec<_>>();

        let _pc_result = Self::prepared_batch_check_evaluations(
            cs,
            prepared_verification_key,
            &lc_commitments,
            &query_set,
            &evaluations,
            proofs,
            opening_challenges,
            opening_challenges_bits,
            batching_rands,
            batching_rands_bits,
        )?;

        Ok(())
    }

    fn create_labeled_commitment_gadget(
        label: String,
        commitment: Self::CommitmentGadget,
        degree_bound: Option<UInt32>,
    ) -> Self::LabeledCommitmentGadget {
        Self::LabeledCommitmentGadget {
            label,
            commitment,
            degree_bound,
        }
    }

    fn create_prepared_labeled_commitment_gadget(
        label: String,
        prepared_commitment: Self::PreparedCommitmentGadget,
        degree_bound: Option<UInt32>,
    ) -> Self::PreparedLabeledCommitmentGadget {
        Self::PreparedLabeledCommitmentGadget {
            label,
            prepared_commitment,
            degree_bound,
            np_phantom: PhantomData,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BTreeMap;
    use crate::{Evaluations, LabeledPolynomial, PolynomialCommitment, QuerySet};
    use algebra::{
        mnt4_298::MNT4_298,
        mnt4_298::{Fq, Fr},
        mnt6_298::MNT6_298,
        test_rng,
    };
    use algebra_core::{One, UniformRand, Zero};
    use nonnative::params::MNT64Small as TestNonNativeFieldParams;
    use r1cs_core::ConstraintSystem;
    use r1cs_std::{
        alloc::AllocGadget, mnt6_298::PairingGadget as MNT6PairingGadget,
        test_constraint_system::TestConstraintSystem,
    };
    use rand::rngs::mock::StepRng;
    use rand::{distributions::Distribution, Rng};

    pub(crate) use ff_fft::DensePolynomial as Polynomial;

    type TestCommitmentScheme = MarlinKZG10<MNT6_298>;
    type TestCheckGadget =
        MarlinKZG10Gadget<MNT298Cycle, TestNonNativeFieldParams, MNT6PairingGadget>;
    type TestCommitmentGadget =
        LabeledCommitmentGadget<MNT298Cycle, TestNonNativeFieldParams, MNT6PairingGadget>;
    type TestProofGadget = ProofGadget<MNT298Cycle, TestNonNativeFieldParams, MNT6PairingGadget>;
    type TestVKGadget = VerifierKeyGadget<MNT298Cycle, TestNonNativeFieldParams, MNT6PairingGadget>;
    type TestNonNativeGadget = NonNativeFieldGadget<TestNonNativeFieldParams>;

    #[derive(Copy, Clone, Debug)]
    struct MNT298Cycle;

    impl CycleEngine for MNT298Cycle {
        type E1 = MNT4_298;
        type E2 = MNT6_298;
    }

    fn marlin_pc_batch_check_test(prepared: bool) {
        let num_inputs = 100;
        //let mut rng = &mut thread_rng();
        let rng = &mut StepRng::new(0, 2);
        let mut inputs: Vec<Option<Fr>> = Vec::with_capacity(num_inputs);
        for _ in 0..num_inputs {
            inputs.push(Some(rng.gen()));
        }

        {
            let mut cs = TestConstraintSystem::<Fr>::new();

            let max_degree = 30;
            let supported_degree = rand::distributions::Uniform::from(1..=max_degree).sample(rng);

            let pp = TestCommitmentScheme::setup(max_degree, rng).unwrap();

            let num_points_in_query_set = 3;
            let num_polynomials = 5;
            let opening_challenge = Fq::rand(rng);

            let mut opening_challenges = Vec::new();
            let mut batching_rands = Vec::new();

            let mut labels = Vec::new();

            let mut degree_bounds = Some(Vec::new());

            let mut polynomials = Vec::new();

            let mut num_of_polys_with_degree_bounds = 0;
            for i in 0..num_polynomials {
                let label = format!("Test{}", i);
                labels.push(label.clone());
                let degree = rand::distributions::Uniform::from(1..=supported_degree).sample(rng);
                let poly = Polynomial::rand(degree, rng);

                let degree_bound = if let Some(degree_bounds) = &mut degree_bounds {
                    if rng.gen() {
                        let range = rand::distributions::Uniform::from(degree..=max_degree);
                        let degree_bound = range.sample(rng);
                        degree_bounds.push(degree_bound);
                        num_of_polys_with_degree_bounds += 1;
                        Some(degree_bound)
                    } else {
                        None
                    }
                } else {
                    None
                };
                let hiding_bound = None;

                polynomials.push(LabeledPolynomial::new_owned(
                    label,
                    poly,
                    degree_bound,
                    hiding_bound,
                ))
            }

            let mut current_challenge = Fq::one();
            for _ in 0..num_polynomials + num_of_polys_with_degree_bounds {
                opening_challenges.push(current_challenge);
                current_challenge *= &opening_challenge;
            }

            for _ in 0..num_points_in_query_set {
                batching_rands.push(Fq::rand(rng));
            }

            let (ck, vk) = TestCommitmentScheme::trim(
                &pp,
                supported_degree,
                degree_bounds.as_ref().map(|s| s.as_slice()),
            )
            .unwrap();

            let (comms, rands) =
                TestCommitmentScheme::commit(&ck, &polynomials, Some(rng)).unwrap();


            let mut opening_challenge_gadgets = Vec::new();
            let mut opening_challenge_bits = Vec::new();
            for (i, challenge) in opening_challenges.iter().enumerate() {
                let challenge_gadget = TestNonNativeGadget::alloc_input(
                    cs.ns(|| format!("opening challenge {}", i)),
                    || Ok(challenge.clone()),
                )
                .unwrap();

                let challenge_bits = challenge_gadget
                    .to_bits(cs.ns(|| format!("opening challenge {} to bits", i)))
                    .unwrap();

                opening_challenge_gadgets.push(challenge_gadget);
                opening_challenge_bits.push(challenge_bits);
            }

            let mut batching_rands_gadgets = Vec::new();
            let mut batching_rands_bits = Vec::new();
            for (i, batching_rand) in batching_rands.iter().enumerate() {
                let batching_rands_gadget = TestNonNativeGadget::alloc_input(
                    cs.ns(|| format!("batching rand {}", i)),
                    || Ok(batching_rand.clone()),
                )
                .unwrap();

                let bits = batching_rands_gadget
                    .to_bits(cs.ns(|| format!("batching rand {} to bits", i)))
                    .unwrap();

                batching_rands_gadgets.push(batching_rands_gadget);
                batching_rands_bits.push(bits);
            }

            // Construct query set
            let mut query_set = QuerySet::new();
            let mut evaluations = Evaluations::new();
            let mut query_set_gadget = QuerySetGadget::<TestNonNativeFieldParams>::new();
            let mut evaluations_gadget = EvaluationsGadget::<TestNonNativeFieldParams>::new();
            for i in 0..num_points_in_query_set {
                let mut point = Fq::one();
                point += point;
                point += point;
                point += Fq::one();
                point += point;
                point += Fq::one();

                let point_gadget =
                    TestNonNativeGadget::alloc(cs.ns(|| format!("point {}", i)), || Ok(point))
                        .unwrap();

                for (j, label) in labels.iter().enumerate() {
                    query_set.insert((label.to_string(), point));
                    query_set_gadget.insert((label.to_string(), point_gadget.clone()));

                    let value = polynomials[j].evaluate(point);
                    let value_gadget = TestNonNativeGadget::alloc_input(
                        cs.ns(|| format!("value {} {}", i, j)),
                        || Ok(value.clone()),
                    )
                    .unwrap();

                    evaluations.insert((label.to_string(), point), value);
                    evaluations_gadget
                        .insert((label.to_string(), point_gadget.clone()), value_gadget);
                }
            }

            let mut commitment_gadgets: Vec<TestCommitmentGadget> = Vec::new();
            for (c, comm) in comms.iter().enumerate() {
                let commitment_gadget =
                    TestCommitmentGadget::alloc(cs.ns(|| format!("commitment {}", c)), || {
                        Ok(comm.clone())
                    })
                    .unwrap();
                commitment_gadgets.push(commitment_gadget);
            }

            let batch_proof = TestCommitmentScheme::batch_open(
                &ck,
                &polynomials,
                &comms,
                &query_set,
                opening_challenge.clone(),
                &rands,
                Some(rng),
            )
            .unwrap();
            let proof_vec: Vec<_> = batch_proof.clone().into();
            let proof_gadgets: Vec<_> = proof_vec
                .iter()
                .enumerate()
                .map(|(i, proof)| {
                    TestProofGadget::alloc(cs.ns(|| format!("proof {}", i)), || {
                        Ok(proof.clone())
                    })
                    .unwrap()
                })
                .collect();

            // Check native proof
            let native_result = TestCommitmentScheme::batch_check(
                &vk,
                &comms,
                &query_set,
                &evaluations,
                &batch_proof,
                opening_challenge.clone(),
                rng,
            )
            .unwrap();
            assert!(native_result);

            let vk_gadget =
                TestVKGadget::alloc_input(cs.ns(|| "verifier key"), || Ok(vk.clone())).unwrap();

            // Check proof in constraints world
            if prepared {
                let prepared_vk_gadget = PreparedVerifierKeyGadget::prepare(
                    &mut cs.ns(|| "prepare vk gadget"),
                    &vk_gadget,
                )
                .unwrap();

                let prepared_commitment_gadgets = commitment_gadgets
                    .iter()
                    .enumerate()
                    .map(|(c, comm)| {
                        PreparedLabeledCommitmentGadget::prepare(
                            &mut cs.ns(|| format!("prepare commitment {}", c)),
                            &comm,
                        )
                        .unwrap()
                    })
                    .collect();

                <TestCheckGadget as PCCheckGadget<
                    Fq,
                    TestCommitmentScheme,
                    TestNonNativeFieldParams,
                    Fr,
                >>::prepared_batch_check_evaluations(
                    cs.ns(|| "checking evaluation"),
                    &prepared_vk_gadget,
                    &prepared_commitment_gadgets,
                    &query_set_gadget,
                    &evaluations_gadget,
                    &proof_gadgets,
                    &opening_challenge_gadgets,
                    &opening_challenge_bits,
                    &batching_rands_gadgets,
                    &batching_rands_bits,
                )
                .unwrap();
            } else {
                <TestCheckGadget as PCCheckGadget<
                    Fq,
                    TestCommitmentScheme,
                    TestNonNativeFieldParams,
                    Fr,
                >>::batch_check_evaluations(
                    cs.ns(|| "checking evaluation"),
                    &vk_gadget,
                    &commitment_gadgets,
                    &query_set_gadget,
                    &evaluations_gadget,
                    &proof_gadgets,
                    &opening_challenge_gadgets,
                    &opening_challenge_bits,
                    &batching_rands_gadgets,
                    &batching_rands_bits,
                )
                .unwrap();
            }

            if !cs.is_satisfied() {
                println!("\n=========================================================");
                println!("\nUnsatisfied constraints:");
                println!("\n{:?}", cs.which_is_unsatisfied().unwrap());
                println!("\n=========================================================");
            }
            assert!(cs.is_satisfied());

            // Constraints profiling
            let show_constraints = false;
            if show_constraints {
                let constraint_str_list = cs.get_constraints_list();
                let prefix_vec = vec![
                    "setup",
                    "verifier key",
                    "checking evaluation/process batch 0",
                    "checking evaluation/process batch 1",
                    "checking evaluation/process batch 2",
                    "checking evaluation/random initializer to subtract",
                    "checking evaluation/batch check/initialize total_c",
                    "checking evaluation/batch check/initialize total_w",
                    "checking evaluation/batch check/process combined commitment 0",
                    "checking evaluation/batch check/process combined commitment 1",
                    "checking evaluation/batch check/process combined commitment 2",
                    "checking evaluation/batch check/prepare for pairing check",
                    "checking evaluation/batch check/lhs = pairing(total_w, beta_h) * pairing(total_c, h)",
                    "checking evaluation/batch check/lhs = rhs",
                ];

                let mut map: BTreeMap<String, (u64, u64, u64, u64)> = BTreeMap::new();

                for name in &constraint_str_list {
                    let mut flag = false;
                    let mut matched_prefix: Option<String> = None;

                    let name = name.replace("#", "/");

                    for prefix in &prefix_vec {
                        let pattern = format!("{}/", prefix);
                        if name.starts_with(pattern.as_str()) {
                            flag = true;
                            matched_prefix = Some(String::from(*prefix));
                            break;
                        }
                    }

                    let key = match flag {
                        true => matched_prefix.unwrap(),
                        false => name.clone(),
                    };

                    if map.contains_key(&key) {
                        let val = map.get(&key).unwrap().clone();

                        if name.contains("debug_marker_counting_number_of_muls") {
                            map.insert(key.clone(), (val.0, val.1, val.2 + 1, val.3));
                        } else if name.contains("debug_marker_counting_number_of_reduces") {
                            map.insert(key.clone(), (val.0, val.1, val.2, val.3 + 1));
                        } else if name.contains("debug_marker_counting_number_of_adds") {
                            map.insert(key.clone(), (val.0, val.1 + 1, val.2, val.3));
                        } else {
                            map.insert(key.clone(), (val.0 + 1, val.1, val.2, val.3));
                        }
                    } else {
                        let val = (0, 0, 0, 0);

                        if name.contains("debug_marker_counting_number_of_muls") {
                            map.insert(key.clone(), (val.0, val.1, val.2 + 1, val.3));
                        } else if name.contains("debug_marker_counting_number_of_reduces") {
                            map.insert(key.clone(), (val.0, val.1, val.2, val.3 + 1));
                        } else if name.contains("debug_marker_counting_number_of_adds") {
                            map.insert(key.clone(), (val.0, val.1 + 1, val.2, val.3));
                        } else {
                            map.insert(key.clone(), (val.0 + 1, val.1, val.2, val.3));
                        }
                    }
                }

                let mut vec: Vec<(String, (u64, u64, u64, u64))> = Vec::with_capacity(map.len());

                for (k, v) in map.iter() {
                    vec.push((k.clone(), v.clone()));
                }

                vec.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

                for (k, v) in vec.iter() {
                    println!(
                        "{} has {} constraints, {} adds, {} muls, and {} reduces",
                        k, v.0, v.1, v.2, v.3
                    );
                }
            }
        }
    }

    #[test]
    fn marlin_pc_batch_check_test_unprepared() {
        marlin_pc_batch_check_test(false)
    }

    #[test]
    fn marlin_pc_batch_check_test_prepared() {
        marlin_pc_batch_check_test(true)
    }

    #[test]
    fn marlin_pc_equation_check_test() {
        let num_iters = 1;
        let max_degree = Some(10);
        let supported_degree = Some(2);
        let num_polynomials = 3;
        let enforce_degree_bounds = true;
        let max_num_queries = 3;
        let num_equations = Some(1);

        let mut cs = TestConstraintSystem::<Fr>::new();

        let rng = &mut test_rng();
        let max_degree =
            max_degree.unwrap_or(rand::distributions::Uniform::from(2..=64).sample(rng));

        let pp = TestCommitmentScheme::setup(max_degree, rng).unwrap();

        for iteration in 0..num_iters {
            let mut cs = cs.ns(|| format!("iteration {}", iteration));

            let supported_degree = supported_degree
                .unwrap_or(rand::distributions::Uniform::from(1..=max_degree).sample(rng));
            assert!(
                max_degree >= supported_degree,
                "max_degree < supported_degree"
            );
            let mut polynomials = Vec::new();
            let mut degree_bounds = if enforce_degree_bounds {
                Some(Vec::new())
            } else {
                None
            };

            let opening_challenge = Fq::rand(rng);

            let mut labels = Vec::new();
            println!("Sampled supported degree");

            // Generate polynomials
            let num_points_in_query_set =
                rand::distributions::Uniform::from(1..=max_num_queries).sample(rng);
            for i in 0..num_polynomials {
                let label = format!("Test{}", i);
                labels.push(label.clone());
                let degree = rand::distributions::Uniform::from(1..=supported_degree).sample(rng);
                let poly = Polynomial::rand(degree, rng);

                let degree_bound = if let Some(degree_bounds) = &mut degree_bounds {
                    if rng.gen() {
                        let range = rand::distributions::Uniform::from(degree..=supported_degree);
                        let degree_bound = range.sample(rng);
                        degree_bounds.push(degree_bound);

                        Some(degree_bound)
                    } else {
                        None
                    }
                } else {
                    None
                };

                let hiding_bound = None;
                println!("Hiding bound: {:?}", hiding_bound);

                polynomials.push(LabeledPolynomial::new_owned(
                    label,
                    poly,
                    degree_bound,
                    hiding_bound,
                ))
            }

            println!("supported degree: {:?}", supported_degree);
            println!("num_points_in_query_set: {:?}", num_points_in_query_set);
            println!("{:?}", degree_bounds);
            println!("{}", num_polynomials);
            println!("{}", enforce_degree_bounds);

            let (ck, vk) = TestCommitmentScheme::trim(
                &pp,
                supported_degree,
                degree_bounds.as_ref().map(|s| s.as_slice()),
            )
            .unwrap();
            println!("Trimmed");

            let (comms, rands) =
                TestCommitmentScheme::commit(&ck, &polynomials, Some(rng)).unwrap();

            // Let's construct our equations
            let mut linear_combinations = Vec::new();
            let mut query_set = QuerySet::new();
            let mut evaluations = Evaluations::new();
            let mut query_set_gadget = QuerySetGadget::<TestNonNativeFieldParams>::new();
            let mut evaluations_gadget = EvaluationsGadget::<TestNonNativeFieldParams>::new();

            let mut opening_challenges = Vec::new();
            let mut batching_rands = Vec::new();
            let mut max_num_opening_challenges_needed = 0;

            for i in 0..num_points_in_query_set {
                let point = Fq::rand(rng);
                let point_gadget =
                    TestNonNativeGadget::alloc(cs.ns(|| format!("point {}", i)), || Ok(point))
                        .unwrap();
                let mut num_opening_challenges_needed = 0;

                for j in 0..num_equations.unwrap() {
                    let label = format!("query {} eqn {}", i, j);
                    let mut lc = LinearCombination::empty(label.clone());

                    let mut value = Fq::zero();
                    let mut has_degree_bound = false;
                    let should_have_degree_bounds: bool = rng.gen();
                    for (k, label) in labels.iter().enumerate() {
                        if should_have_degree_bounds {
                            value += &polynomials[k].evaluate(point.clone());
                            lc.push((Fq::one(), label.to_string().into()));

                            if polynomials[k].degree_bound().is_some() {
                                has_degree_bound = true;
                            }

                            break;
                        } else {
                            let poly = &polynomials[k];
                            if poly.degree_bound().is_some() {
                                continue;
                            } else {
                                assert!(poly.degree_bound().is_none());
                                let coeff = Fq::rand(rng);
                                value += &(coeff * poly.evaluate(point.clone()));
                                lc.push((coeff.clone(), label.to_string().into()));
                            }
                        }
                    }

                    num_opening_challenges_needed += 1;
                    if has_degree_bound == true {
                        num_opening_challenges_needed += 1;
                    }

                    let value_gadget = TestNonNativeGadget::alloc_input(
                        cs.ns(|| format!("value {} {}", i, j)),
                        || Ok(value.clone()),
                    )
                    .unwrap();

                    evaluations.insert((label.clone(), point.clone()), value);
                    evaluations_gadget
                        .insert((label.to_string(), point_gadget.clone()), value_gadget);

                    if !lc.is_empty() {
                        linear_combinations.push(lc);
                        // Insert query
                        query_set.insert((label.clone(), point.clone()));
                        query_set_gadget.insert((label.clone(), point_gadget.clone()));
                    }
                }

                if num_opening_challenges_needed > max_num_opening_challenges_needed {
                    max_num_opening_challenges_needed = num_opening_challenges_needed;
                }
            }

            let mut current_challenge = Fq::one();
            for _ in 0..max_num_opening_challenges_needed {
                opening_challenges.push(current_challenge);
                current_challenge *= &opening_challenge;
            }

            for _ in 0..num_points_in_query_set {
                batching_rands.push(Fq::rand(rng));
            }

            if linear_combinations.is_empty() {
                continue;
            }
            println!("Generated query set");
            println!("Linear combinations: {:?}", linear_combinations);

            let mut commitment_gadgets: Vec<TestCommitmentGadget> = Vec::new();
            for (c, comm) in comms.iter().enumerate() {
                let commitment_gadget =
                    TestCommitmentGadget::alloc(cs.ns(|| format!("commitment {}", c)), || {
                        Ok(comm.clone())
                    })
                    .unwrap();
                commitment_gadgets.push(commitment_gadget);
            }

            let mut opening_challenge_gadgets = Vec::new();
            let mut opening_challenge_bits = Vec::new();
            for (i, challenge) in opening_challenges.iter().enumerate() {
                let challenge_gadget = TestNonNativeGadget::alloc_input(
                    cs.ns(|| format!("opening challenge {}", i)),
                    || Ok(challenge.clone()),
                )
                .unwrap();

                let challenge_bits = challenge_gadget
                    .to_bits(cs.ns(|| format!("opening challenge {} to bits", i)))
                    .unwrap();

                opening_challenge_gadgets.push(challenge_gadget);
                opening_challenge_bits.push(challenge_bits);
            }

            let mut batching_rands_gadgets = Vec::new();
            let mut batching_rands_bits = Vec::new();
            for (i, batching_rand) in batching_rands.iter().enumerate() {
                let batching_rand_gadget = TestNonNativeGadget::alloc_input(
                    cs.ns(|| format!("batching rand {}", i)),
                    || Ok(batching_rand.clone()),
                )
                .unwrap();

                let batching_rand_bits = batching_rand_gadget
                    .to_bits(cs.ns(|| format!("batching rand {} to bits", i)))
                    .unwrap();

                batching_rands_gadgets.push(batching_rand_gadget);
                batching_rands_bits.push(batching_rand_bits);
            }

            let proof = TestCommitmentScheme::open_combinations(
                &ck,
                &linear_combinations,
                &polynomials,
                &comms,
                &query_set,
                opening_challenge.clone(),
                &rands,
                Some(rng),
            )
            .unwrap();
            println!("Generated proof");

            // Check native proof
            let native_result = TestCommitmentScheme::check_combinations(
                &vk,
                &linear_combinations,
                &comms,
                &query_set,
                &evaluations,
                &proof,
                opening_challenge.clone(),
                rng,
            )
            .unwrap();
            assert!(native_result);

            // Convert outputs from native world to constraints world

            // Convert proof
            let proof_vec: Vec<_> = proof.proof.clone().into();
            let proof_gadgets: Vec<_> = proof_vec
                .iter()
                .enumerate()
                .map(|(i, proof)| {
                    TestProofGadget::alloc(cs.ns(|| format!("proof {}", i)), || Ok(proof.clone()))
                        .unwrap()
                })
                .collect();
            let evals_vec: Vec<_> = query_set_gadget
                .iter()
                .map(|q| evaluations_gadget.get(q).unwrap().clone())
                .collect();
            let batch_lc_proof_gadget = BatchLCProofGadget {
                proofs: proof_gadgets,
                evals: Some(evals_vec),
            };

            // Convert linear combinations
            let linear_combination_gadgets = linear_combinations
                .iter()
                .enumerate()
                .map(|(l, lc)| {
                    LinearCombinationGadget::alloc(
                        cs.ns(|| format!("linear combination {}", l)),
                        || Ok(lc),
                    )
                    .unwrap()
                })
                .collect();

            // Convert verifier key
            let vk_gadget =
                TestVKGadget::alloc_input(cs.ns(|| "verifier key"), || Ok(vk.clone())).unwrap();
            let prepared_vk_gadget =
                PreparedVerifierKeyGadget::prepare(&mut cs.ns(|| "prepare vk gadget"), &vk_gadget)
                    .unwrap();

            // Prepare commitment gadgets
            let prepared_commitment_gadgets = commitment_gadgets
                .iter()
                .enumerate()
                .map(|(c, comm)| {
                    PreparedLabeledCommitmentGadget::prepare(
                        &mut cs.ns(|| format!("prepare commitment {}", c)),
                        &comm,
                    )
                    .unwrap()
                })
                .collect();

            // Check proof in constraints world
            <TestCheckGadget as PCCheckGadget<
                Fq,
                TestCommitmentScheme,
                TestNonNativeFieldParams,
                Fr,
            >>::prepared_check_combinations(
                cs.ns(|| "checking evaluation"),
                &prepared_vk_gadget,
                &linear_combination_gadgets,
                &prepared_commitment_gadgets,
                &query_set_gadget,
                &evaluations_gadget,
                &batch_lc_proof_gadget,
                &opening_challenge_gadgets,
                &opening_challenge_bits,
                &batching_rands_gadgets,
                &batching_rands_bits,
            )
            .unwrap();
        }

        if !cs.is_satisfied() {
            println!("\n=========================================================");
            println!("\nUnsatisfied constraints:");
            println!("\n{:?}", cs.which_is_unsatisfied().unwrap());
            println!("\n=========================================================");
        }
        assert!(cs.is_satisfied());
    }
}
