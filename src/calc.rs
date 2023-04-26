extern crate nalgebra as na;

use crate::Dtype;
use na::*;

/// Linear equations: A*x = b,
/// In this case: A for static_kmat, x for disps, b for forces
pub struct LinearEqs<const N: usize> {
    state: bool,
    pub disps: [Dtype; N],
    pub forces: [Dtype; N],
    pub static_kmat: [[Dtype; N]; N],
}

impl<const N: usize> LinearEqs<N> {
    pub fn new(disps: [Dtype; N], forces: [Dtype; N], static_kmat: [[Dtype; N]; N]) -> Self {
        LinearEqs {
            disps,
            forces,
            state: false,
            static_kmat,
        }
    }

    /// get disp on every single node
    pub fn disps_rlt(&mut self) -> &[Dtype; N] {
        if self.state == false {
            self.lu_direct_solver();
            &self.disps
        } else {
            &self.disps
        }
    }

    /// get force on every single node
    pub fn forces_rlt(&mut self) -> &[Dtype; N] {
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
        let kmat = SMatrix::<Dtype, N, N>::from(self.static_kmat).transpose();
        let force = SVector::from(self.forces);

        let disps_unknown_idx = nonzero_disps_idx(&self.disps);
        let force_known = force.select_rows(disps_unknown_idx.iter());
        let kmat_eff = kmat
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());

        // solve the K.q = F by LU decomposition
        let disps_unknown: Vec<Dtype> = kmat_eff.lu().solve(&force_known).unwrap().data.into();
        println!("\n>>> LU decomposition method down!\n\tresult: err = 0\n");

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
    pub fn gauss_seidel_iter_solver(&mut self, calc_error: Dtype) {
        // pre-process
        let unknown_disps_idx = nonzero_disps_idx(&self.disps);

        let force = SVector::from(self.forces);
        let f_eff = force.select_rows(unknown_disps_idx.iter());

        let kmat = SMatrix::<Dtype, N, N>::from(self.static_kmat).transpose();
        let kmat_eff = kmat
            .select_columns(unknown_disps_idx.iter())
            .select_rows(unknown_disps_idx.iter());

        // construct Gauss-Seidel iter method
        let l = kmat_eff.lower_triangle();
        let d = DMatrix::from_diagonal(&kmat_eff.diagonal());
        let u = kmat_eff.upper_triangle() - &d;
        let grad = -l.try_inverse().unwrap();

        let size = unknown_disps_idx.len();
        let mut x = DVector::<Dtype>::zeros(size);

        // Gauss-Seidel iterator loop
        let mut count: usize = 0;
        loop {
            let tmp = &grad * &u * &x - &grad * &f_eff;
            //println!("#{}, x={}, err={}", &count, &x, (&tmp - &x).abs().max());

            if (&tmp - &x).abs().max() < calc_error {
                println!(
                    "\n>>> Gauss-Seidel iter method down!\n\tresult: iter = {},\n\t\terr  = {:8.6}\n",
                    count,
                    (&tmp - &x).abs().max()
                );
                break;
            }
            x = tmp;
            count += 1usize;
        }

        // write result
        let _: Vec<_> = unknown_disps_idx
            .iter()
            .enumerate()
            .map(|(i, &idx)| self.disps[idx] = x[i])
            .collect();
        self.forces = (((kmat * SVector::from(self.disps)) - force) + force).into();
        self.state = true;
    }
}

/// 将Kmat中节点位移已知的自由度找出来
fn nonzero_disps_idx<'a, T: IntoIterator<Item = &'a Dtype>>(container: T) -> Vec<usize> {
    let idx: Vec<usize> = container
        .into_iter()
        .enumerate()
        .filter(|(_, &ele)| ele != 0.0)
        .map(|(idx, _)| idx)
        .collect();
    idx
}
