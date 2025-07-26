use crate::dtty::basic::Dtype;

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
