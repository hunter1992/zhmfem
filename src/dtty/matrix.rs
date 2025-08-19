use crate::dtty::basic::Dtype;
use std::boxed::Box;

/// Symmetric Skyline method to storage element and part stiffness matrix
#[derive(Clone)]
pub struct CompressedMatrix {
    pub values: Box<Vec<Dtype>>,
    pub pointr: Box<Vec<usize>>,
}

impl CompressedMatrix {
    /// Construct a compressed matrix
    pub fn new(values: Box<Vec<Dtype>>, pointr: Box<Vec<usize>>) -> Self {
        CompressedMatrix { values, pointr }
    }

    /// Recover a square matrix from compressed mat which is in CompressedMatrix shape
    pub fn recover<const DIM: usize>(&self) -> Box<[[Dtype; DIM]; DIM]> {
        let mut matrix: Box<[[Dtype; DIM]; DIM]> = Box::new([[0.0; DIM]; DIM]);
        for idx in 0..DIM {
            let len = self.pointr[idx + 1] - self.pointr[idx];
            for jump in 0..len {
                matrix[idx][idx - jump] = self.values[self.pointr[idx] + len - jump - 1];
                matrix[idx - jump][idx] = self.values[self.pointr[idx] + len - jump - 1];
            }
        }
        matrix
    }
}
