/// This file define all kinds of data structure of ZHMFEM
use atomic_float::{AtomicF32, AtomicF64};
use na::SMatrix;

pub type Dtype = f32;
pub type ADtype = AtomicF32;

pub type Jacobian2D = SMatrix<Dtype, 2, 2>;
pub type Jacobian3D = SMatrix<Dtype, 3, 3>;

/// Store strain or stress data
#[derive(Clone)]
pub enum Data {
    Dim1(Vec<[Dtype; 1]>),
    Dim2(Vec<[Dtype; 3]>),
    Dim3(Vec<[Dtype; 6]>),
}

/// Symmetric Skyline method to storage element and part stiffness matrix
#[derive(Clone)]
pub struct CompressedMatrix {
    pub value: Vec<Dtype>,
    pub ptr: Vec<usize>,
}

impl CompressedMatrix {
    pub fn new(value: Vec<Dtype>, ptr: Vec<usize>) -> Self {
        CompressedMatrix { value, ptr }
    }

    pub fn recover<const D: usize>(&self) -> [[Dtype; D]; D] {
        let mut arr: [[Dtype; D]; D] = [[0.0; D]; D];
        for idx in 0..D {
            let len = self.ptr[idx + 1] - self.ptr[idx];
            for jump in 0..len {
                arr[idx][idx - jump] = self.value[self.ptr[idx] + len - jump - 1];
                arr[idx - jump][idx] = self.value[self.ptr[idx] + len - jump - 1];
            }
        }
        arr
    }
}
