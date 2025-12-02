use crate::dtty::{
    basic::Dtype,
    matrix::{CompressedMatrixCSR, CompressedMatrixSKS},
};
use crate::tool::compress_matrix_csr_1based;
use na::{DMatrix, DVector, SMatrix, SVector};
use std::collections::HashSet;
use std::time::Instant;

pub struct LinearEqs<const D: usize> {
    state: bool,
    pub disps: [Dtype; D],
    pub loads: [Dtype; D],
    pub disps_0_idx: Vec<usize>,
    pub static_kmat: CompressedMatrixSKS,
    pub external_force: Option<[Dtype; D]>,
    pub solver_time_consuming: Option<std::time::Duration>,
}

impl<const D: usize> LinearEqs<D> {
    pub fn new(
        disps: [Dtype; D],
        loads: [Dtype; D],
        disps_0_idx: Vec<usize>,
        static_kmat: CompressedMatrixSKS,
    ) -> Self {
        LinearEqs {
            state: false,
            disps,
            loads,
            disps_0_idx,
            static_kmat,
            external_force: None,
            solver_time_consuming: None,
        }
    }

    /// Solve the defined problem
    pub fn solve(&mut self, solve_method: &str, calc_error: Dtype) {
        match solve_method {
            "lu" => self.lu_direct_solver(),
            // "pardiso" => self.panua_pardiso_direct_solver(),
            "cholesky" => self.cholesky_direct_solver(),
            "gs" => self.gauss_seidel_iter_solver(calc_error),
            _ => panic!("!!! The provided algorithm is currently not supported!"),
        }
    }

    /// 使用Panua Tech的PARDISO进行求解
    // pub fn panua_pardiso_direct_solver(&mut self) {
    //     // pre-process
    //     let kmat: CompressedMatrixCSR = compress_matrix_csr(
    //         SMatrix::<Dtype, D, D>::from(*self.static_kmat.recover()).transpose(),
    //     );
    //     let loads = SVector::from(self.loads);

    //     let disps_unknown_idx = idx_subtract::<D>(self.disps_0_idx.clone());
    //     let force_known = loads.select_rows(disps_unknown_idx.iter());
    //     let kmat_eff = kmat
    //         .select_columns(disps_unknown_idx.iter())
    //         .select_rows(disps_unknown_idx.iter());
    // }

    /// 使用LU分解求解线性方程组(直接法)
    /// solve "A * x = b" using LU decomposition method
    pub fn lu_direct_solver(&mut self) {
        // pre-process
        let kmat = SMatrix::<Dtype, D, D>::from(*self.static_kmat.recover()).transpose();
        let loads = SVector::from(self.loads);

        let disps_unknown_idx = idx_subtract::<D>(self.disps_0_idx.clone());
        let force_known = loads.select_rows(disps_unknown_idx.iter());
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
        self.external_force = Some((kmat * SVector::from(self.disps)).into());
        self.solver_time_consuming = Some(duration_lu);
        self.state = true;
    }

    pub fn cholesky_direct_solver(&mut self) {
        // pre-process
        let kmat = SMatrix::<Dtype, D, D>::from(*self.static_kmat.recover()).transpose();
        let loads = SVector::from(self.loads);

        let disps_unknown_idx = idx_subtract::<D>(self.disps_0_idx.clone());
        let force_known = loads.select_rows(disps_unknown_idx.iter());
        let kmat_eff = kmat
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());

        let time_lu = Instant::now();
        // solve the K.q = F by LU decomposition
        let disps_unknown_rlt: Vec<Dtype> =
            kmat_eff.cholesky().unwrap().solve(&force_known).data.into();
        let duration_lu = time_lu.elapsed();
        println!(
            "\n>>> Cholesky decomposition method down!\n\ttime consuming = {:?}",
            duration_lu
        );

        // write result into fields
        let _: Vec<_> = disps_unknown_idx
            .iter()
            .enumerate()
            .map(|(i, &idx)| self.disps[idx] = disps_unknown_rlt[i])
            .collect();
        self.external_force = Some((kmat * SVector::from(self.disps)).into());
        self.solver_time_consuming = Some(duration_lu);
        self.state = true;
    }

    /// 使用Guass-Seidel迭代求解线性方程组(数值方法)
    /// solve "A*x = F" using Gauss-Seidel iterator:
    ///   x(k+1) = -[(D+L)^(-1)] * U * x(k)  + [(D+L)^(-1)] * F
    pub fn gauss_seidel_iter_solver(&mut self, calc_error: Dtype) {
        // pre-process
        let disps_unknown_idx = idx_subtract::<D>(self.disps_0_idx.clone());

        let loads = SVector::from(self.loads);
        let f_eff = loads.select_rows(disps_unknown_idx.iter());

        let kmat = SMatrix::<Dtype, D, D>::from(*self.static_kmat.recover()).transpose();
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
                self.solver_time_consuming = Some(duration_gs);
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
        self.external_force = Some((kmat * SVector::from(self.disps)).into());
        self.state = true;
    }
}

/// 求出Part各单元所有节点位移向量中，非零位移的索引
fn idx_subtract<const N: usize>(zero_disps_idx: Vec<usize>) -> Vec<usize> {
    let all_idx: [usize; N] = std::array::from_fn(|x| x);
    let whole: HashSet<usize> = HashSet::from(all_idx);
    let zeros: HashSet<usize> = zero_disps_idx.into_iter().collect();
    (&whole - &zeros).iter().cloned().collect()
}
