use crate::dtty::basic::Dtype;
use std::fmt;

/// Symmetric Skyline method to storage element and part stiffness matrix
#[derive(Clone)]
pub struct CompressedMatrixSKS {
    pub values: Vec<Dtype>, // matrix element data
    pub pointr: Vec<usize>, //
}

impl CompressedMatrixSKS {
    /// Construct a compressed matrix
    pub fn new(values: Vec<Dtype>, pointr: Vec<usize>) -> Self {
        CompressedMatrixSKS { values, pointr }
    }

    /// Recover a square matrix from compressed mat which is in CompressedMatrix shape
    pub fn recover<const DIM: usize>(&self) -> Box<[[Dtype; DIM]; DIM]> {
        let mut matrix: Box<[[Dtype; DIM]; DIM]> = Box::new([[0.0; DIM]; DIM]);
        for row_idx in 0..DIM {
            let len = self.pointr[row_idx + 1] - self.pointr[row_idx];
            for jump in 0..len {
                matrix[row_idx][row_idx - jump] =
                    self.values[self.pointr[row_idx] + len - jump - 1];
                matrix[row_idx - jump][row_idx] =
                    self.values[self.pointr[row_idx] + len - jump - 1];
            }
        }
        matrix
    }
}

impl fmt::Display for CompressedMatrixSKS {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\n\tdata: {:?}\n\tpointer: {:?}\n",
            self.values, self.pointr
        )
    }
}

/// Compressed sparse row format to storage element's or part's stiffness matrix
/// Notice: here the stiffness matrix is symmetric, so the compressed matrix in csr
/// format just hava half data.
#[derive(Clone)]
pub struct CompressedMatrixCSR {
    pub values: Vec<Dtype>, // matrix element data
    pub colidx: Vec<usize>, // data in which colum
    pub pointr: Vec<usize>, // rows head's pointer
}

impl CompressedMatrixCSR {
    /// Construct a compressed matrix
    pub fn new(values: Vec<Dtype>, colidx: Vec<usize>, pointr: Vec<usize>) -> Self {
        CompressedMatrixCSR {
            values,
            colidx,
            pointr,
        }
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

impl fmt::Display for CompressedMatrixCSR {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\n\tdata: {:?}\n\tcolumn index: {:?}\n\tpointer: {:?}\n",
            self.values, self.colidx, self.pointr
        )
    }
}
