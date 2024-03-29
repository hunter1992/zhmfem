extern crate nalgebra as na;

use crate::Dtype;
use na::*;
use std::collections::HashSet;
use std::time::Instant;

/// Non-Linear equations: A(x) * x = b，这里的A(x)是Kt，也就是切线刚度
/// 矩阵，是未知的节点位移向量的函数，右侧的b是节点受到的外力，是保守力。

/// Linear equations: A*x = b,
/// In this case: A for known static_kmat, x for disps, b for forces
pub struct LinearEqs<const N: usize> {
    pub disps: [Dtype; N],
    pub forces: [Dtype; N],
    pub disps_0_idx: Vec<usize>,
    pub static_kmat: [[Dtype; N]; N],
    state: bool, // solved or not flag
}

impl<const N: usize> LinearEqs<N> {
    pub fn new(
        disps: [Dtype; N],
        forces: [Dtype; N],
        disps_0_idx: Vec<usize>,
        static_kmat: [[Dtype; N]; N],
    ) -> Self {
        LinearEqs {
            disps,
            forces,
            disps_0_idx,
            static_kmat,
            state: false,
        }
    }

    /// get disp on every single node
    pub fn disps_rlt(&mut self, solve_method: usize, calc_error: Dtype) -> &[Dtype; N] {
        if self.state == false {
            if solve_method == 0 {
                self.lu_direct_solver();
            } else if solve_method == 1 {
                self.gauss_seidel_iter_solver(calc_error);
            }
            &self.disps
        } else {
            &self.disps
        }
    }

    /// get force on every single node
    pub fn forces_rlt(&mut self, solve_method: usize, calc_error: Dtype) -> &[Dtype; N] {
        if self.state == false {
            if solve_method == 0 {
                self.lu_direct_solver();
            } else if solve_method == 1 {
                self.gauss_seidel_iter_solver(calc_error);
            }
            &self.forces
        } else {
            &self.forces
        }
    }

    /// 使用LU分解求解线性方程组(直接法)
    /// solve "A * x = b" using LU decomposition method
    pub fn lu_direct_solver(&mut self) {
        // pre-process
        let kmat = SMatrix::<Dtype, N, N>::from(self.static_kmat).transpose();
        let force = SVector::from(self.forces);

        //let disps_unknown_idx = nonzero_disps_idx(&self.disps); //旧算法停用
        let disps_unknown_idx = idx_subtract::<N>(self.disps_0_idx.clone());
        let force_known = force.select_rows(disps_unknown_idx.iter());
        let kmat_eff = kmat
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());

        let time_lu = Instant::now();
        // solve the K.q = F by LU decomposition
        let disps_unknown_rlt: Vec<Dtype> = kmat_eff.lu().solve(&force_known).unwrap().data.into();
        let duration_lu = time_lu.elapsed();
        println!(
            "\n>>> LU decomposition method down!\n\ttime consuming = {:?}",
            duration_lu
        );

        // write result into fields
        let _: Vec<_> = disps_unknown_idx
            .iter()
            .enumerate()
            .map(|(i, &idx)| self.disps[idx] = disps_unknown_rlt[i])
            .collect();
        self.forces = (((kmat * SVector::from(self.disps)) - force) + force).into();
        self.state = true;
    }

    /// 使用Guass-Seidel迭代求解线性方程组(数值方法)
    /// solve "A*x = F" using Gauss-Seidel iterator:
    ///   x(k+1) = -[(D+L)^(-1)] * U * x(k)  + [(D+L)^(-1)] * F
    pub fn gauss_seidel_iter_solver(&mut self, calc_error: Dtype) {
        // pre-process
        //let disps_unknown_idx = nonzero_disps_idx(&self.disps); //旧算法停用
        let disps_unknown_idx = idx_subtract::<N>(self.disps_0_idx.clone());

        let force = SVector::from(self.forces);
        let f_eff = force.select_rows(disps_unknown_idx.iter());

        let kmat = SMatrix::<Dtype, N, N>::from(self.static_kmat).transpose();
        let kmat_eff = kmat
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());

        // construct Gauss-Seidel iter method
        let l = kmat_eff.lower_triangle();
        let d = DMatrix::from_diagonal(&kmat_eff.diagonal());
        let u = kmat_eff.upper_triangle() - &d;
        let grad = -l.try_inverse().unwrap();

        let size = disps_unknown_idx.len();
        let mut x = DVector::<Dtype>::zeros(size);

        // Gauss-Seidel iterator loop
        let mut count: usize = 0;
        let time_gs = Instant::now();
        loop {
            let tmp = &grad * &u * &x - &grad * &f_eff;
            //println!("#{}, x={}, err={}", &count, &x, (&tmp - &x).abs().max());

            if (&tmp - &x).abs().max() < calc_error {
                let duration_gs = time_gs.elapsed();
                print!("\n>>> Gauss-Seidel iter method down!");
                println!(
                    "\n\ttime consuming = {:?}\n\tresult:   iter = {}\n\t\t   err = {:8.6}",
                    duration_gs,
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
fn nonzero_disps_idx<'a, T: IntoIterator<Item = &'a Dtype>>(container: T) -> Vec<usize> {
    let idx: Vec<usize> = container
        .into_iter()
        .enumerate()
        .filter(|(_, &ele)| ele != 0.0)
        .map(|(idx, _)| idx)
        .collect();
    idx
}
/// 求出Part中所有节点位移向量中，非零位移的索引
fn idx_subtract<const N: usize>(zero_disps_idx: Vec<usize>) -> Vec<usize> {
    let all_idx: [usize; N] = std::array::from_fn(|x| x);
    let whole: HashSet<usize> = HashSet::from(all_idx);
    let zeros: HashSet<usize> = zero_disps_idx.into_iter().collect();
    (&whole - &zeros).iter().cloned().collect()
}
