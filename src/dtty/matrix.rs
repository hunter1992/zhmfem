use crate::dtty::{aligndata::AlignedVec, basic::Dtype};
use crate::tool::{compress_symmetry_matrix_csr, compress_symmetry_matrix_sks};
use std::fmt;

/// Symmetric Skyline method to storage element and part stiffness matrix
/// For example, there is a sparse symmetric matrix:
/// (comes from Panua Tech's PARDISO user guide v8.2 page8)
///
/// [[ 7,   0,   1,   0,   0,   2,   7,   0]
///  [ 0,  -4,   8,   0,   2,   0,   0,   0]
///  [ 1,   8,   1,   0,   0,   0,   0,   5]
///  [ 0,   0,   0,   7,   0,   0,   9,   0]
///  [ 0,   2,   0,   0,   5,  -1,   5,   0]
///  [ 2,   0,   0,   0,  -1,   0,   0,   5]
///  [ 7,   0,   0,   9,   5,   0,  11,   0]
///  [ 0,   0,   5,   0,   0,   5,   0,   5]]
///
/// The values storage:
/// [7  -4   1  8  1  7  2  0  0  5  2  0  0  0  -1  0  7  0  0  9  5  0  11  5  0  0  5  0  5]
///  ^   ^   ^        ^  ^           ^                  ^                     ^
/// [0   1   2        5  6           10                 16                    23  29]
/// Previous line of code is Pointers to columns.
#[derive(Clone)]
pub struct CompressedMatrixSKS {
    pub values: Vec<Dtype>, // matrix element data
    pub pointr: Vec<usize>, //
}

impl CompressedMatrixSKS {
    /// Compresse a matrix in SKS form from a full square matrix array
    pub fn new<const DIM: usize>(square_mat: &[[Dtype; DIM]; DIM]) -> Self {
        compress_symmetry_matrix_sks(square_mat)
    }

    /// Compresse a matrix in SKS form from a full square matrix Vec
    pub fn from_vec(mat: Vec<Vec<Dtype>>) -> Self {
        assert_eq!(mat.len(), mat[0].len());
        let dim = mat.len();
        let mut values: Vec<Dtype> = vec![];
        let mut pointr: Vec<usize> = vec![];

        let mut flag: bool = true; //当前行是否需要处理
        let mut col_idx_in_single_row: usize = 0;
        for row_idx in 0..dim {
            pointr.push(col_idx_in_single_row);
            for col_idx in 0..=row_idx {
                if mat[row_idx][col_idx] == 0.0 {
                    if flag == false {
                        // 正定矩阵必满秩，下面处理全零行(或列)的分支暂时注释掉
                        /*if idy == idx {
                            value.push(mat[idx][idy]);
                            counter += 1;
                        }*/
                        continue;
                    } else {
                        values.push(mat[row_idx][col_idx]);
                        col_idx_in_single_row += 1;
                        continue;
                    }
                }
                values.push(mat[row_idx][col_idx]);
                col_idx_in_single_row += 1;
                flag = true;
            }
            flag = false;
        }
        pointr.push(values.len());

        CompressedMatrixSKS { values, pointr }
    }

    /// Recover a square matrix from compressed mat which is in CompressedMatrix shape
    pub fn recover_square_arr<const DIM: usize>(&self) -> [[Dtype; DIM]; DIM] {
        let mut matrix: [[Dtype; DIM]; DIM] = [[0.0; DIM]; DIM];
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

    pub fn recover_square_vec<const DIM: usize>(&self) -> Vec<Vec<Dtype>> {
        let mut matrix: Vec<Vec<Dtype>> = vec![vec![0.0; DIM]; DIM];
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

    /// Get the element at (x, y) in original square matrix
    #[inline]
    pub fn get_elem(&self, idx: usize, idy: usize) -> Dtype {
        if idx < idy {
            // 压缩的时候存的下三角
            return self.get_elem(idy, idx);
        }
        let band_length: usize = self.pointr[idx + 1] - self.pointr[idx];
        let empty_width: usize = idx + 1 - band_length;
        if idy < empty_width {
            0.0
        } else {
            self.values[self.pointr[idx] - empty_width + idy]
        }
    }

    /// Generate a sub-matrix from CompressedMatrixSKS
    pub fn get_sub_matrix(&self, index: &[usize]) -> Vec<Vec<Dtype>> {
        let dim: usize = index.len();
        let mut sub_mat: Vec<Vec<Dtype>> = vec![vec![0.0; dim]; dim];
        for idx in 0..dim {
            for idy in 0..dim {
                sub_mat[idx][idy] = self.get_elem(index[idx], index[idy]);
            }
        }
        sub_mat
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
///
/// [[ 7,   0,   1,   0,   0,   2,   7,   0]
///  [ 0,  -4,   8,   0,   2,   0,   0,   0]
///  [ 1,   8,   1,   0,   0,   0,   0,   5]
///  [ 0,   0,   0,   7,   0,   0,   9,   0]
///  [ 0,   2,   0,   0,   5,  -1,   5,   0]
///  [ 2,   0,   0,   0,  -1,   0,   0,   5]
///  [ 7,   0,   0,   9,   5,   0,  11,   0]
///  [ 0,   0,   5,   0,   0,   5,   0,   5]]
///
///
/// data:         [7.0, 1.0, 2.0, 7.0, -4.0, 8.0, 2.0, 1.0, 5.0, 7.0, 9.0, 5.0, -1.0, 5.0, 5.0, 11.0, 5.0]
/// column index: [0,   2,   5,   6,    1,   2,   4,   2,   7,   3,   6,   4,    5,   6,   7,   6,    7  ]
/// row pointer:  [0,                   4,             7,        9,        11,             14,  15,   16, 17]

#[derive(Clone)]
pub struct CompressedMatrixCSR {
    pub values: Vec<Dtype>, // matrix element data
    pub colidx: Vec<usize>, // data in which colum
    pub rowptr: Vec<usize>, // rows head's pointer
    pub baseinfo: u8,
    pub squaredim: usize,
}

impl CompressedMatrixCSR {
    /// Compresse a matrix in CSR form from a full square matrix array
    pub fn new<const DIM: usize>(square_mat: &[[Dtype; DIM]; DIM]) -> Self {
        compress_symmetry_matrix_csr(square_mat)
    }

    /// Compresse a matrix in CSR form from a full square matrix Vec
    pub fn from_vec(mat: Vec<Vec<Dtype>>) -> Self {
        let mut values: Vec<Dtype> = vec![];
        let mut colidx: Vec<usize> = vec![];
        let mut rowptr: Vec<usize> = vec![];
        let baseinfo: u8 = 0;
        let squaredim: usize = mat.len();

        let mut no_head_in_current_row: bool = true;
        let mut counter: usize = 0;
        for row_idx in 0..squaredim {
            for col_idx in row_idx..squaredim {
                //右上三角
                if 0.0 == mat[row_idx][col_idx] {
                    continue;
                } else {
                    values.push(mat[row_idx][col_idx]);
                    colidx.push(col_idx);
                    if no_head_in_current_row {
                        rowptr.push(counter);
                        no_head_in_current_row = false;
                    }
                    counter += 1;
                }
            }
            no_head_in_current_row = true;
        }

        rowptr.push(counter);

        CompressedMatrixCSR {
            values,
            colidx,
            rowptr,
            baseinfo,
            squaredim,
        }
    }

    /// Convert 0-based CSR matrix to 1-based
    pub fn convert_1base(&mut self) {
        if self.baseinfo == 0 {
            for idx in self.colidx.iter_mut() {
                *idx += 1;
            }
            for ptr in self.rowptr.iter_mut() {
                *ptr += 1;
            }
            self.baseinfo = 1;
        } else if self.baseinfo == 1 {
            panic!(
                "!!! The CSR matrix is already 1-based, from src/dtty/matrix/CompressedMatrixCSR.convert_1base."
            )
        }
    }

    /// Convert 1-based CSR matrix to 0-based
    pub fn convert_0base(&mut self) {
        if self.baseinfo == 1 {
            self.colidx
                .iter_mut()
                .for_each(|column_idx| *column_idx -= 1);
            self.rowptr
                .iter_mut()
                .for_each(|rowpointer| *rowpointer -= 1);
            self.baseinfo = 0;
        } else if self.baseinfo == 0 {
            panic!(
                "!!! The CSR matrix is already 0-based, from src/dtty/matrix/CompressedMatrixCSR.convert_0base."
            )
        }
    }

    /// Recover a square matrix from compressed mat which is in CompressedMatrix shape
    pub fn recover_square_arr<const DIM: usize>(&mut self) -> Box<[[Dtype; DIM]; DIM]> {
        let mut matrix: Box<[[Dtype; DIM]; DIM]> = Box::new([[0.0; DIM]; DIM]);
        if self.baseinfo == 1 {
            self.convert_0base();
        }
        if self.colidx.len() != *self.rowptr.last().unwrap() {
            panic!("!!! Error from CompressedMatrixCSR.recover, wrong length.")
        }
        for current_row in 0..DIM {
            let n_elems_in_current_row = self.rowptr[current_row + 1] - self.rowptr[current_row];
            for n in 0..n_elems_in_current_row {
                matrix[current_row][self.colidx[self.rowptr[current_row] + n]] =
                    self.values[self.rowptr[current_row] + n];
                matrix[self.colidx[self.rowptr[current_row] + n]][current_row] =
                    self.values[self.rowptr[current_row] + n];
            }
        }
        matrix
    }

    pub fn recover_square_vec<const DIM: usize>(&mut self) -> Vec<Vec<Dtype>> {
        let mut matrix: Vec<Vec<Dtype>> = vec![vec![0.0; DIM]; DIM];
        if self.baseinfo == 1 {
            self.convert_0base();
        }
        if self.colidx.len() != *self.rowptr.last().unwrap() {
            panic!("!!! Error from CompressedMatrixCSR.recover, wrong length.")
        }
        for current_row in 0..DIM {
            let n_elems_in_current_row = self.rowptr[current_row + 1] - self.rowptr[current_row];
            for n in 0..n_elems_in_current_row {
                matrix[current_row][self.colidx[self.rowptr[current_row] + n]] =
                    self.values[self.rowptr[current_row] + n];
                matrix[self.colidx[self.rowptr[current_row] + n]][current_row] =
                    self.values[self.rowptr[current_row] + n];
            }
        }
        matrix
    }

    /// Generate a CSR format matrix that meets Panua PARDISO requirements
    pub fn as_panua_pardiso_arg(self) -> CSRforPanuaPARDISO {
        let colidx: Vec<i32> = self.colidx.into_iter().map(|val| val as i32).collect();
        let rowptr: Vec<i32> = self.rowptr.into_iter().map(|val| val as i32).collect();
        CSRforPanuaPARDISO {
            values: AlignedVec::<f64>::from_slice(&self.values),
            colidx: AlignedVec::<i32>::from_slice(&colidx),
            rowptr: AlignedVec::<i32>::from_slice(&rowptr),
            squaredim: self.squaredim,
        }
    }
}

impl fmt::Display for CompressedMatrixCSR {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\n\tdata: {:?}\n\tcolumn index: {:?}\n\tpointer: {:?}\n",
            self.values, self.colidx, self.rowptr
        )
    }
}

/// The CSR matrix used to serve as the PARDISO parameter of
/// Panua Tech is slightly different from the usual CSR matrix:
/// 1) CSR for Panua's PARDISO is 1-based;
/// 2) colum index & row pointer element is i32 not usize
pub struct CSRforPanuaPARDISO {
    pub values: AlignedVec<f64>, // matrix element data
    pub colidx: AlignedVec<i32>, // data in which colum
    pub rowptr: AlignedVec<i32>, // rows head's pointer
    pub squaredim: usize,
}

impl CSRforPanuaPARDISO {
    /// Construct from 0_based CompressedMatrixCSR
    pub fn new<const DIM: usize>(square_mat: &[[Dtype; DIM]; DIM]) -> Self {
        let mut csr: CompressedMatrixCSR = compress_symmetry_matrix_csr(square_mat);
        csr.convert_1base();
        csr.as_panua_pardiso_arg()
    }

    /// 对 CSR 矩阵的所有组件进行内存预热
    #[inline(always)]
    pub fn warm_up(&self) {
        #[cfg(feature = "parallel")]
        {
            // 只有开启了特性，且当前 CPU 核心数 > 1 时才使用并行
            // 这样可以兼顾编译速度和运行时的极端性能
            rayon::join(
                || self.warm_up_slice(&self.values),
                || {
                    rayon::join(
                        || self.warm_up_slice(&self.colidx),
                        || self.warm_up_slice(&self.rowptr),
                    )
                },
            );
        }

        #[cfg(not(feature = "parallel"))]
        {
            // 默认走超轻量级的单线程预热，编译极快
            self.warm_up_slice(&self.values);
            self.warm_up_slice(&self.colidx);
            self.warm_up_slice(&self.rowptr);
        }
    }

    /// 内部辅助函数：以 64 字节（Cache Line）为步长进行强制读取
    #[inline(always)]
    pub fn warm_up_slice<T>(&self, data: &[T]) {
        if data.is_empty() {
            return;
        }

        let ptr = data.as_ptr() as *const u8;
        let size = data.len() * std::mem::size_of::<T>();
        let mut offset = 0;

        unsafe {
            while offset < size {
                // 使用 volatile 读取每个缓存行的第一个字节
                // 这会强制触发 MMU 填充 TLB，并让硬件预取器感知到流式访问
                std::ptr::read_volatile(ptr.add(offset));
                // 现代 CPU 缓存行通常为 64 字节，步进 64 可实现最高效率覆盖
                offset += 64;
            }

            // 确保最后一个字节也被触碰到，防止大对象末尾页缺失
            if size > 0 {
                std::ptr::read_volatile(ptr.add(size - 1));
            }
        }
    }

    /// Compresse a matrix in CSR form from a full square matrix Vec
    pub fn from_vec(mat: Vec<Vec<Dtype>>) -> Self {
        let mut values: Vec<Dtype> = vec![];
        let mut colidx: Vec<i32> = vec![];
        let mut rowptr: Vec<i32> = vec![];
        let squaredim: usize = mat.len();

        let mut no_head_in_current_row: bool = true;
        let mut counter: usize = 0;
        for row_idx in 0..squaredim {
            for col_idx in row_idx..squaredim {
                if 0.0 == mat[row_idx][col_idx] {
                    continue;
                } else {
                    values.push(mat[row_idx][col_idx]);
                    colidx.push(col_idx as i32);
                    if no_head_in_current_row {
                        rowptr.push(counter as i32);
                        no_head_in_current_row = false;
                    }
                    counter += 1;
                }
            }
            no_head_in_current_row = true;
        }

        rowptr.push(counter as i32);

        for idx in colidx.iter_mut() {
            *idx += 1;
        }
        for ptr in rowptr.iter_mut() {
            *ptr += 1;
        }

        CSRforPanuaPARDISO {
            values: AlignedVec::<f64>::from_slice(&values),
            colidx: AlignedVec::<i32>::from_slice(&colidx),
            rowptr: AlignedVec::<i32>::from_slice(&rowptr),
            squaredim,
        }
    }

    pub fn to_normal_csr(self) -> CompressedMatrixCSR {
        let mut csr = CompressedMatrixCSR {
            values: self.values.into_vec(),
            colidx: self.colidx.into_vec(),
            rowptr: self.rowptr.into_vec(),
            baseinfo: 1,
            squaredim: self.squaredim,
        };
        csr.convert_0base();
        csr
    }
}

impl fmt::Display for CSRforPanuaPARDISO {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\n\tdata: {:?}\n\tcolumn index: {:?}\n\tpointer: {:?}\n\tbase_from: {}\n",
            self.values, self.colidx, self.rowptr, self.squaredim
        )
    }
}
