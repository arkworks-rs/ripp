//#![deny(warnings, unused, future_incompatible, nonstandard_style)]

use std::{
    error::Error as ErrorTrait,
    fmt::{Display, Formatter, Result as FmtResult},
};

// pub mod applications;
pub mod gipa;
pub mod ip_commitment;
// pub mod tipa;

pub type Error = Box<dyn ErrorTrait>;

#[derive(Debug)]
pub enum InnerProductArgumentError {
    MessageLengthInvalid(usize, usize),
    InnerProductInvalid,
}

impl ErrorTrait for InnerProductArgumentError {
    fn source(self: &Self) -> Option<&(dyn ErrorTrait + 'static)> {
        None
    }
}

impl Display for InnerProductArgumentError {
    fn fmt(self: &Self, f: &mut Formatter<'_>) -> FmtResult {
        let msg = match self {
            InnerProductArgumentError::MessageLengthInvalid(left, right) => {
                format!("left length, right length: {left}, {right}")
            }
            InnerProductArgumentError::InnerProductInvalid => "inner product not sound".to_string(),
        };
        write!(f, "{}", msg)
    }
}
