use crate::dtty::{
    basic::Dtype,
    matrix::{CSRforPanuaPARDISO, CompressedMatrixCSR, CompressedMatrixSKS},
};
use crate::tool::idx_subtract;
use na::{DMatrix, DVector, SMatrix, SVector};
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

    /// Solving systems of linear algebraic equations
    /// 1) solve_method:
    ///                  auto     ---- Panua Tech's PARDISO
    ///                  lu       ---- LU decomposition method
    ///                  pardiso  ---- Panua Tech's PARDISO
    ///                  cholesky ---- Cholesky decomposition method
    ///                  gs       ---- Gauss-Seidel iterative algorithm
    ///
    /// 2) calc_error  : set the truncation error for the iterative algorithm
    ///
    /// 3) cpus        : set the number of threads for parallel computing
    pub fn solve(
        &mut self,
        solve_algorithm: &str,
        truncation_error: Dtype,
        cpus: usize,
    ) -> Result<bool, i32> {
        match solve_algorithm {
            "pardiso" => self.panua_pardiso_direct_solver(cpus),
            "cholesky" => self.cholesky_direct_solver(),
            "lu" => self.lu_direct_solver(),
            "gs" => self.gauss_seidel_iter_solver(truncation_error),
            "auto" => {
                if 0 < D && D <= 256 {
                    self.gauss_seidel_iter_solver(truncation_error)
                } else if 256 < D && D <= 2048 {
                    self.cholesky_direct_solver()
                } else {
                    self.panua_pardiso_direct_solver(cpus)
                }
            }
            _ => panic!("No solution algorithm selected! If unsure, you can enter \"auto\"."),
        }
    }

    /// Solve the "A*x = b" with direct method using Panua Tech's PARDISO
    pub fn panua_pardiso_direct_solver(&mut self, cpus: usize) -> Result<bool, i32> {
        use crate::calc::panua;
        use std::env;
        use std::os::raw::c_void;
        use std::ptr;

        // STEP 0. Prepare the coefficient matrix and RHS vector needed for calculation
        let disps_unknown_idx: Vec<usize> = idx_subtract::<D>(self.disps_0_idx.clone());
        let mut b: Vec<Dtype> = disps_unknown_idx
            .iter()
            .map(|&idx| self.loads[idx])
            .collect();
        let mut csr: CompressedMatrixCSR =
            CompressedMatrixCSR::from_vec(self.static_kmat.get_sub_matrix(&disps_unknown_idx));

        csr.convert_1base();
        let matrix: CSRforPanuaPARDISO = csr.as_panua_pardiso_arg();

        let mut x: Vec<Dtype> = vec![0.0; matrix.squaredim];

        // STEP 1. Set environment variables and shared library file paths
        unsafe {
            env::set_var("PARDISO_LIC_PATH", "/opt/panua-pardiso-20240229-linux/lic");
            env::set_var("PARDISO_PATH", "/opt/panua-pardiso-20240229-linux/lib");
            env::set_var("PARDISOLICMESSAGE", "1"); // 抑制许可证检查信息输出
        }
        let pardiso_path = env::var("PARDISO_PATH").expect("!!! PARDISO_PATH must be set!");
        let pardiso_lic_path =
            env::var("PARDISO_LIC_PATH").expect("!!! PARDISO_LIC_PATH must be set!");
        println!("\n>>> PARDISO_PATH = {:?}", pardiso_path);
        println!("\n>>> PARDISO_LIC_PATH= {:?}", pardiso_lic_path);

        // STEP 2. Set PARDISO parameters
        const NRHS: usize = 1;
        let mut pt: [*mut c_void; 64] = [ptr::null_mut(); 64];
        let mut maxfct: i32 = 1;
        let mut mnum: i32 = 1;
        let mut mtype: i32 = 2;
        let mut phase: i32;
        let mut n: i32 = matrix.squaredim as i32;
        let a: Vec<Dtype> = matrix.values;
        let ia: Vec<i32> = matrix.rowptr;
        let ja: Vec<i32> = matrix.colidx;
        let mut nrhs: i32 = NRHS as i32;
        let mut iparm: [i32; 64] = [0; 64];
        let mut msglvl: i32 = 1;
        let mut dparm: [f64; 64] = [0.0; 64];
        let mut error: i32 = 0;
        let mut solver: i32 = 0;

        unsafe {
            // ---------------------
            // Initialization
            // ---------------------
            iparm[0] = 0;
            panua::pardisoinit(
                pt.as_mut_ptr(),
                &mut mtype,
                &mut solver,
                iparm.as_mut_ptr(),
                dparm.as_mut_ptr(),
                &mut error,
            );
            if error != 0 {
                println!("\n[PARDISO]: !!!ERROR during Reordering");
                return Err(error);
            } else {
                println!("\n[PARDISO]: Initialization complete.");
            }

            // ---------------------
            // Ordering
            // ---------------------
            iparm[2] = cpus as i32;
            phase = 11;
            let time_start1 = Instant::now();
            panua::pardiso(
                pt.as_mut_ptr(),
                &mut maxfct,
                &mut mnum,
                &mut mtype,
                &mut phase,
                &mut n,
                a.as_ptr() as *mut f64,
                ia.as_ptr() as *mut i32,
                ja.as_ptr() as *mut i32,
                ptr::null_mut(),
                &mut nrhs,
                iparm.as_mut_ptr(),
                &mut msglvl,
                ptr::null_mut(),
                ptr::null_mut(),
                &mut error,
                dparm.as_mut_ptr(),
            );
            let time_end1 = time_start1.elapsed();
            if error != 0 {
                println!("\n[PARDISO]: !!!ERROR during Reordering");
                return Err(error);
            }

            // ------------------------
            // Numerical Factorization
            // ------------------------
            phase = 22;
            let time_start2 = Instant::now();
            panua::pardiso(
                pt.as_mut_ptr(),
                &mut maxfct,
                &mut mnum,
                &mut mtype,
                &mut phase,
                &mut n,
                a.as_ptr() as *mut f64,
                ia.as_ptr() as *mut i32,
                ja.as_ptr() as *mut i32,
                ptr::null_mut(),
                &mut nrhs,
                iparm.as_mut_ptr(),
                &mut msglvl,
                ptr::null_mut(),
                ptr::null_mut(),
                &mut error,
                dparm.as_mut_ptr(),
            );
            let time_end2 = time_start2.elapsed();
            if error != 0 {
                return Err(error);
            }

            // ------------------------
            // Solving
            // ------------------------
            phase = 33;
            let time_start3 = Instant::now();
            panua::pardiso(
                pt.as_mut_ptr(),
                &mut maxfct,
                &mut mnum,
                &mut mtype,
                &mut phase,
                &mut n,
                a.as_ptr() as *mut f64,
                ia.as_ptr() as *mut i32,
                ja.as_ptr() as *mut i32,
                ptr::null_mut(),
                &mut nrhs,
                iparm.as_mut_ptr(),
                &mut msglvl,
                b.as_mut_ptr(),
                x.as_mut_ptr(),
                &mut error,
                dparm.as_mut_ptr(),
            );
            let time_end3 = time_start3.elapsed();
            let time_end4 = time_start1.elapsed();
            if error != 0 {
                return Err(error);
            }

            println!("\n>>> PARDISO calculation down!");
            println!(
                "\tN of equations:                       {}",
                disps_unknown_idx.len()
            );
            println!("\tAnalysis time consuming:              {:?}", time_end1);
            println!("\tNumeric Factorization time consuming: {:?}", time_end2);
            println!("\tSolve time consuming:                 {:?}", time_end3);
            println!("\tTotal time consuming:                 {:?}", time_end4);
            self.solver_time_consuming = Some(time_end4);
        }

        // write result into fields
        let _: Vec<_> = disps_unknown_idx
            .iter()
            .enumerate()
            .map(|(i, &idx)| self.disps[idx] = x[i])
            .collect();
        let static_kmat = SMatrix::<Dtype, D, D>::from(self.static_kmat.recover_square_arr());
        self.external_force = Some((static_kmat * SVector::from(self.disps)).into());
        self.state = true;

        Ok(self.state)
    }

    /// Solve "A * x = b" with direct method using LU decomposition
    pub fn lu_direct_solver(&mut self) -> Result<bool, i32> {
        // pre-process
        let static_kmat = SMatrix::<Dtype, D, D>::from(self.static_kmat.recover_square_arr());
        let loads = SVector::from(self.loads);

        let disps_unknown_idx = idx_subtract::<D>(self.disps_0_idx.clone());
        let force_known = loads.select_rows(disps_unknown_idx.iter());
        let kmat_eff = static_kmat
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());

        let time_start_lu = Instant::now();
        // solve the K.q = F by LU decomposition
        let disps_unknown_rlt: Vec<Dtype> = kmat_eff.lu().solve(&force_known).unwrap().data.into();
        let duration_lu = time_start_lu.elapsed();
        println!(
            "\n>>> LU decomposition method down!\n\tN of equations:  {}\n\tTime consuming: {:?}",
            disps_unknown_idx.len(),
            duration_lu
        );

        // write result into fields
        let _: Vec<_> = disps_unknown_idx
            .iter()
            .enumerate()
            .map(|(i, &idx)| self.disps[idx] = disps_unknown_rlt[i])
            .collect();
        self.external_force = Some((static_kmat * SVector::from(self.disps)).into());
        self.solver_time_consuming = Some(duration_lu);
        self.state = true;

        Ok(self.state)
    }

    /// Solve "A * x = b" with direct method using Cholesky decomposition
    pub fn cholesky_direct_solver(&mut self) -> Result<bool, i32> {
        // pre-process
        let static_kmat = SMatrix::<Dtype, D, D>::from(self.static_kmat.recover_square_arr());
        let loads = SVector::from(self.loads);

        let disps_unknown_idx = idx_subtract::<D>(self.disps_0_idx.clone());
        let force_known = loads.select_rows(disps_unknown_idx.iter());
        let kmat_eff = static_kmat
            .select_columns(disps_unknown_idx.iter())
            .select_rows(disps_unknown_idx.iter());

        let time_start_cholesky = Instant::now();
        // solve the K.q = F by LU decomposition
        let disps_unknown_rlt: Vec<Dtype> =
            kmat_eff.cholesky().unwrap().solve(&force_known).data.into();
        let duration_lu = time_start_cholesky.elapsed();
        println!(
            "\n>>> Cholesky decomposition method down!\n\tN of equations: {}\n\tTime consuming: {:?}",
            disps_unknown_idx.len(),
            duration_lu
        );

        // write result into fields
        let _: Vec<_> = disps_unknown_idx
            .iter()
            .enumerate()
            .map(|(i, &idx)| self.disps[idx] = disps_unknown_rlt[i])
            .collect();
        self.external_force = Some((static_kmat * SVector::from(self.disps)).into());
        self.solver_time_consuming = Some(duration_lu);
        self.state = true;

        Ok(self.state)
    }

    /// Solve "A*x = F" with iterator method using Gauss-Seidel
    ///   x(k+1) = -[(D+L)^(-1)] * U * x(k)  + [(D+L)^(-1)] * F
    /// url:https://www.jianshu.com/p/e14d9e910984
    pub fn gauss_seidel_iter_solver(&mut self, calc_error: Dtype) -> Result<bool, i32> {
        // pre-process
        let disps_unknown_idx: Vec<usize> = idx_subtract::<D>(self.disps_0_idx.clone());

        let loads = SVector::from(self.loads);
        let f_eff = loads.select_rows(disps_unknown_idx.iter());

        let static_kmat = SMatrix::<Dtype, D, D>::from(self.static_kmat.recover_square_arr());
        let kmat_eff = static_kmat
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
        self.external_force = Some((static_kmat * SVector::from(self.disps)).into());
        self.state = true;

        Ok(self.state)
    }
}
