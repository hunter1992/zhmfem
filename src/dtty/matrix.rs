use crate::dtty::basic::Dtype;
use crate::tool::compress_matrix_sks;
use std::fmt;

/// Symmetric Skyline method to storage element and part stiffness matrix
/// For example, there is a sparse symmetric matrix:
///
/// [[ 11,  12,   0,  14,   0,   0,   0,   0,   0]
///  [      22,  23,  24,  25,   0,   0,   0,   0]
///  [           33,  34,   0,   0,  37,   0,   0]
///  [                44,  45,   0,  47,   0,   0]
///  [                     55,  56,  57,   0,  59]
///  [                          66,  67,  68,  69]
///  [                               77,  78,  79]
///  [                                    88,  89]
///  [                                         99]]
///
/// The values storage:
/// [11, 12, 22, 23, 33, 14, 24, 34, 44, 25, 0, 45, 55, 56, 66, 37, 47, 57, 67, 77, 68, 78, 88, 59, 69, 79, 89, 99]
///  ^   ^       ^       ^               ^              ^       ^                   ^           ^
/// [1,  2,      3,      4,              5,             6,      7,                  8,          9]
/// Previous line of code is Pointers to columns.
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

// impl From<CompressedMatrixCSR> for CompressedMatrixSKS {
//     fn from(mat_in_csr: CompressedMatrixCSR) -> Self {
//         let dim: usize = mat_in_csr.pointr.last();
//         compress_matrix_sks(&(mat_in_csr.recover::<dim>()))
//     }
// }

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
        if self.colidx.len() + 1 == self.pointr.len() {
            panic!("!!! Error from CompressedMatrixCSR.recover, wrong length.")
        }
        for idx in 0..DIM {
            let head = self.pointr[idx];
            let tail = self.pointr[idx + 1];
            for i in 0..(head - tail) {
                matrix[idx][self.colidx[i]] = self.values[i];
                matrix[self.colidx[i]][idx] = self.values[i];
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
