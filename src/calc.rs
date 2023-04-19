extern crate nalgebra as na;

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
            self.lu_direct_solver();
            &self.disps
        } else {
            &self.disps
        }
    }

    /// get force on every single node
    pub fn forces_rlt(&mut self) -> &[f64; N] {
        if self.state == false {
            self.lu_direct_solver();
            &self.forces
        } else {
            &self.forces
        }
    }

    /// 使用LU分解求解线性方程组(直接法)  A * x = b
    pub fn lu_direct_solver(&mut self) {
        // pre-process
        let kmat = SMatrix::<f64, N, N>::from(self.static_kmat).transpose();
        let force = SVector::from(self.forces);

        let disps_unknown_idx = nonzero_disps_idx(&self.disps);
        let force_known = force.select_rows(disps_unknown_idx.iter());
        let kmat_eff = kmat
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());

        // solve the K.q = F by LU decomposition
        let disps_unknown: Vec<f64> = kmat_eff.lu().solve(&force_known).unwrap().data.into();

        // write result into fields
        let _: Vec<_> = disps_unknown_idx
            .iter()
            .enumerate()
            .map(|(i, &idx)| self.disps[idx] = disps_unknown[i])
            .collect();
        self.forces = (((kmat * SVector::from(self.disps)) - force) + force).into();
        self.state = true;
    }

    /// 使用Guass-Seidel迭代求解线性方程组(数值方法)
    /// solve 'A*x = F' using Gauss-Seidel iterator:
    ///   x(k+1) = -[(D+L)^(-1)] * U * x(k)  + [(D+L)^(-1)] * F
    pub fn gauss_seidel_iter_solver(&mut self, calc_error: f64) {
        // pre-process
        let force = SVector::from(self.forces);
        let kmat = SMatrix::<f64, N, N>::from(self.static_kmat).transpose();
        let kmat_diag = SMatrix::<f64, N, N>::from(triangle_partition(&self.static_kmat, 0));
        let kmat_up = SMatrix::<f64, N, N>::from(triangle_partition(&self.static_kmat, 1));
        let kmat_low = SMatrix::<f64, N, N>::from(triangle_partition(&self.static_kmat, -1));

        let disps_unknown_idx = nonzero_disps_idx(&self.disps);
        let size = disps_unknown_idx.len();

        // construct Gauss-Seidel iter method
        let f = force.select_rows(disps_unknown_idx.iter());
        let d = kmat_diag
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());
        let l = kmat_low
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());
        let u = kmat_up
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());
        let grad = -&((&d + &l).try_inverse().unwrap());
        let mut x = DVector::<f64>::zeros(size);

        let mut count: usize = 0;
        loop {
            // x(k+1) = -[(D+L)^(-1)] * U * x(k)  + [(D+L)^(-1)] * F
            let tmp = &grad * &u * &x - &grad * &f;
            // println!("#{}, x={}, err={}", &count, &x, (&tmp - &x).abs().max());

            if (&tmp - &x).abs().max() < calc_error {
                println!(
                    ">>> Gauss-Seidel method down!\n\tresult: iter = {},\n\t\terr  = {:8.6}\n",
                    count,
                    (&tmp - &x).abs().max()
                );
                break;
            }
            x = tmp;
            count += 1usize;
        }

        // write result
        let _: Vec<_> = disps_unknown_idx
            .iter()
            .enumerate()
            .map(|(i, &idx)| self.disps[idx] = x[i])
            .collect();
        self.forces = (((kmat * SVector::from(self.disps)) - force) + force).into();
        self.state = true;
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

/// 取方阵的偏移上、下三角阵，或仅含对角元素的方阵
fn triangle_partition<const N: usize>(square_mat: &[[f64; N]; N], mode: isize) -> [[f64; N]; N] {
    let mut rlt: [[f64; N]; N] = [[0.0; N]; N];
    if mode > 0 {
        for i in mode as usize..N {
            for j in 0..i {
                rlt[i][j] = square_mat[i][j];
            }
        }
    } else if mode < 0 {
        for i in (-mode) as usize..N {
            for j in 0..i {
                rlt[j][i] = square_mat[j][i];
            }
        }
    } else {
        for i in 0..N {
            rlt[i][i] = square_mat[i][i];
        }
    }
    rlt
}
