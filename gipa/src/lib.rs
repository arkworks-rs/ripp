use std::{
    error::Error as ErrorTrait,
};

use inner_products::InnerProduct;
use dh_commitments::DoublyHomomorphicCommitment;

pub type Error = Box<dyn ErrorTrait>;


//pub struct GIPA<IP, M1C, M2C, IPC> where
//IP: InnerProduct,
//M1C: DoublyHomomorphicCommitment,
//M2C: DoublyHomomorphicCommitment,
//IPC: DoublyHomomorphicCommitment,
//{
//
//}

