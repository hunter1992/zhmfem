use na::SMatrix;

// Use F32 as basic data type
use atomic_float::AtomicF32;
pub type ADtype = AtomicF32;
pub type Dtype = f32;

// Use F64 as basic data type
// use atomic_float::AtomicF64;
// pub type ADtype = AtomicF64;
// pub type Dtype = f64;

pub type Jacobian2D = SMatrix<Dtype, 2, 2>;
pub type Jacobian3D = SMatrix<Dtype, 3, 3>;
