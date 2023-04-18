extern crate nalgebra as na;

use crate::print_2darr;
use na::*;

/// Linear equations: A*x = b,
/// In this case: A for static_kmat, x for disps, b for forces
pub struct LinearEqs<const N: usize> {
    state: bool,
    pub disps: [f64; N],
    pub forces: [f64; N],
    pub static_kmat: [[f64; N]; N],
}

impl<const N: usize> LinearEqs<N> {
    pub fn new(disps: [f64; N], forces: [f64; N], static_kmat: [[f64; N]; N]) -> Self {
        LinearEqs {
            disps,
            forces,
            state: false,
            static_kmat,
        }
    }

    /// get disp on every single node
    pub fn disps_rlt(&mut self) -> &[f64; N] {
        if self.state == false {
            self.lu_solver();
            &self.disps
        } else {
            &self.disps
        }
    }

    /// get force on every single node
    pub fn forces_rlt(&mut self) -> &[f64; N] {
        if self.state == false {
            self.lu_solver();
            &self.forces
        } else {
            &self.forces
        }
    }

    /// 使用LU分解求解线性方程组(直接法)  A * x = b
    pub fn lu_solver(&mut self) {
        if self.disps.len() != N || self.forces.len() != N {
            panic!("---> Error! from calc::LinearEqs::lu_solver func.");
        }

        let force = SVector::from(self.forces);
        let kmat = SMatrix::<f64, N, N>::from(self.static_kmat);

        let disps_unknown_idx = nonzero_disps_idx(&self.disps);
        let force_known = force.select_rows(disps_unknown_idx.iter());
        let kmat_eff = kmat
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());

        // solve the K.q = F by LU decomposition
        let disps_unknown: Vec<f64> = kmat_eff.lu().solve(&force_known).unwrap().data.into();

        // 写入计算得到的位移和节点力
        let _: Vec<_> = disps_unknown_idx
            .iter()
            .enumerate()
            .map(|(i, &idx)| self.disps[idx] = disps_unknown[i])
            .collect();
        self.forces = (((kmat * SVector::from(self.disps)) - force) + force).into();
        self.state = true;
    }

    /// 使用雅可比迭代求解线性方程组(数值方法)
    pub fn jacobian_iter_solver(&mut self, calc_error: f64) {
        let force = SVector::from(self.forces);
        let kmat = SMatrix::<f64, N, N>::from(self.static_kmat);
        let kmat_diag = SMatrix::<f64, N, N>::from(diag_mat(&self.static_kmat));
        let tri_l = SMatrix::<f64, N, N>::from(tri_l(&self.static_kmat));
        let tri_u = SMatrix::<f64, N, N>::from(tri_u(&self.static_kmat));

        let disps_unknown_idx = nonzero_disps_idx(&self.disps);
        let size = disps_unknown_idx.len();

        /* solve 'A*x = F' using jacobian iterator:
         *     x(k+1) = x(k) * [E - (D^-1)*A] + (D^-1)*F  */

        let E = DMatrix::<f64>::identity(size, size);
        let F = force.select_rows(disps_unknown_idx.iter());
        let A = kmat
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());
        let D = kmat_diag
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());
        let L = tri_l
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());
        let U = tri_u
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());
        let ttmp = (&D + &L).try_inverse().unwrap();
        let B = -&ttmp * U;
        let mut x = DVector::<f64>::zeros(size);

        let mut count: usize = 0;
        loop {
            let tmp = &B * &x + &ttmp * &F;
            println!("#{}, x={}, err={}", &count, &x, (&tmp - &x).abs().max());

            if (&tmp - &x).abs().max() < calc_error {
                break;
            }
            x = tmp;
            count += 1usize;
        }
    }
}

/// 将Kmat中节点位移已知的自由度找出来
fn nonzero_disps_idx<'a, T: IntoIterator<Item = &'a f64>>(container: T) -> Vec<usize> {
    let idx: Vec<usize> = container
        .into_iter()
        .enumerate()
        .filter(|(_, &ele)| ele != 0.0)
        .map(|(idx, _)| idx)
        .collect();
    idx
}

/// 求方阵对角线元素其余元素为0的同尺寸方阵
fn diag_mat<const N: usize>(square_mat: &[[f64; N]; N]) -> [[f64; N]; N] {
    let mut rlt: [[f64; N]; N] = [[0.0; N]; N];
    for i in 0..N {
        for j in 0..N {
            if i == j {
                rlt[i][j] = square_mat[i][j];
            }
        }
    }
    rlt
}
fn tri_l<const N: usize>(square_mat: &[[f64; N]; N]) -> [[f64; N]; N] {
    let mut rlt: [[f64; N]; N] = [[0.0; N]; N];
    for i in 0..N {
        for j in 0..N {
            if i > j {
                rlt[i][j] = square_mat[i][j];
            }
        }
    }
    rlt
}
fn tri_u<const N: usize>(square_mat: &[[f64; N]; N]) -> [[f64; N]; N] {
    let mut rlt: [[f64; N]; N] = [[0.0; N]; N];
    for i in 0..N {
        for j in 0..N {
            if i < j {
                rlt[i][j] = square_mat[i][j];
            }
        }
    }
    rlt
}
