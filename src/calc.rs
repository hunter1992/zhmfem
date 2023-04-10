extern crate nalgebra as na;

use na::*;

/// 线性方程组求解的直接法有LU分解、高斯消去、列主元消去等
/// 但这些方法的缺点在于：对系数矩阵的要求太高，通用性受限。
pub struct Solver<const N: usize> {
    state: bool,
    pub disps: [f64; N],
    pub forces: [f64; N],
    pub static_kmat: [[f64; N]; N],
}

impl<const N: usize> Solver<N> {
    pub fn new(disps: [f64; N], force: [f64; N], static_kmat: [[f64; N]; N]) -> Self {
        Solver {
            disps,
            forces: force,
            state: false,
            static_kmat,
        }
    }

    /// get disp on every single node
    pub fn disps_rlt(&mut self) -> &[f64; N] {
        if self.state == false {
            self.solve_static_lu();
            &self.disps
        } else {
            &self.disps
        }
    }

    /// get force on every single node
    pub fn forces_rlt(&mut self) -> &[f64; N] {
        if self.state == false {
            self.solve_static_lu();
            &self.forces
        } else {
            &self.forces
        }
    }

    /// 使用LU分解求解线性方程组
    pub fn solve_static_lu(&mut self) {
        if self.disps.len() != N || self.forces.len() != N {
            panic!("---> Error! from Solve::solve_static func.");
        }

        let force = SVector::from(self.forces);
        let kmat = SMatrix::<f64, N, N>::from(self.static_kmat);

        let disp_eff_idx = nonzero_disps_idx(&self.disps);
        let force_eff = force.select_rows(disp_eff_idx.iter());
        let kmat_eff = kmat
            .select_columns(disp_eff_idx.iter())
            .select_rows(disp_eff_idx.iter());

        // solve the K.q = F
        let disps_unknown: Vec<f64> = kmat_eff.lu().solve(&force_eff).unwrap().data.into();

        // 写入计算得到的位移和节点力
        let _: Vec<_> = disp_eff_idx
            .iter()
            .enumerate()
            .map(|(i, &idx)| self.disps[idx] = disps_unknown[i])
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
