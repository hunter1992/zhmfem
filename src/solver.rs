extern crate nalgebra as na;

use na::*;

pub struct Solver<const N: usize> {
    status: bool,
    pub disps: [f64; N],
    pub force: [f64; N],
    pub static_kmat: Option<[[f64; N]; N]>,
}

impl<const N: usize> Solver<N> {
    pub fn new(disps: [f64; N], force: [f64; N]) -> Self {
        Solver {
            disps,
            force,
            status: false,
            static_kmat: None,
        }
    }

    /// get disp on every single node
    pub fn get_disps(&mut self) -> &[f64; N] {
        if self.status == false {
            self.solve_static();
            &self.disps
        } else {
            &self.disps
        }
    }

    /// get force on every single node
    pub fn get_forces(&mut self) -> &[f64; N] {
        if self.status == false {
            self.solve_static();
            &self.force
        } else {
            &self.force
        }
    }

    pub fn solve_static(&mut self) {
        if self.disps.len() != N || self.force.len() != N {
            panic!("---> Error! from Solve::solve_static func.");
        }

        let force = SVector::from(self.force);
        let kmat = SMatrix::<f64, N, N>::from(self.static_kmat.unwrap());

        let disp_eff_idx = nonzero_disps_idx(&self.disps);
        let force_eff = force.select_rows(disp_eff_idx.iter());
        let kmat_eff = kmat
            .select_columns(disp_eff_idx.iter())
            .select_rows(disp_eff_idx.iter());

        // solve the K.q = F
        let disps_unknown: Vec<f64> = kmat_eff.lu().solve(&force_eff).unwrap().data.into();
        let _: Vec<_> = disp_eff_idx
            .iter()
            .enumerate()
            .map(|(i, &idx)| self.disps[idx] = disps_unknown[i])
            .collect();
        self.force = (((kmat * SVector::from(self.disps)) - force) + force).into();
    }
}

fn nonzero_disps_idx<'a, T: IntoIterator<Item = &'a f64>>(container: T) -> Vec<usize> {
    let idx: Vec<usize> = container
        .into_iter()
        .enumerate()
        .filter(|(_, &ele)| ele != 0.0)
        .map(|(idx, _)| idx)
        .collect();
    idx
}
