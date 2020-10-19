use crate::Error;
use ark_std::boxed::Box;

#[derive(Debug)]
pub enum LHPCError {
    LHError(Box<dyn ark_std::error::Error>),
    PCError(Box<Error>),
}

impl LHPCError {
    pub fn pc_error(pc_error: Error) -> Self {
        Self::PCError(Box::new(pc_error))
    }

    pub fn lh_error<LHE: 'static + ark_std::error::Error>(lh_error: LHE) -> Self {
        Self::LHError(Box::new(lh_error))
    }
}

impl core::fmt::Display for LHPCError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            LHPCError::LHError(lh_error) => write!(f, "LH Error: {}", (*lh_error)),
            LHPCError::PCError(pc_error) => write!(f, "PC Error: {}", pc_error),
        }
    }
}

impl From<Error> for LHPCError {
    fn from(pc_error: Error) -> Self {
        Self::PCError(Box::new(pc_error))
    }
}

impl ark_std::error::Error for LHPCError {}
