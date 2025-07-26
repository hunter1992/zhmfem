use atomic_float::{AtomicF32, AtomicF64};
use na::SMatrix;

pub type Dtype = f32;
pub type ADtype = AtomicF32;

pub type Jacobian2D = SMatrix<Dtype, 2, 2>;
pub type Jacobian3D = SMatrix<Dtype, 3, 3>;
