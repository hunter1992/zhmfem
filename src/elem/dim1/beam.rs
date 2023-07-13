use crate::{node::Node2D, Dtype, K};
use na::*;
use std::fmt::{self, Write};

/// 细长梁以及纯弯梁，抽象成具有如下特征的理论模型：
///   1) 仅用x坐标刻画所有位置(即细长之含义)；
///   2) 变形主要为垂直于x轴的挠度
/// 但是，节点的自由度不为1，而包含：挠度(y)、转角theta，
/// 此亦是梁单元为非协调单元的原因。
pub struct Beam1D2N<'beam1d2n> {
    pub id: usize,
    pub sec_area: Dtype,
    pub nodes: [&'beam1d2n Node2D; 2],
    pub k_matrix: Option<[[Dtype; 4]; 4]>,
}

impl<'beam1d2n> Beam1D2N<'beam1d2n> {
    /// Generate a Beam1D2N element
    pub fn new(id: usize, sec_area: Dtype, nodes: [&Node2D; 2]) -> Beam1D2N {
        Beam1D2N {
            id,
            sec_area,
            nodes,
            k_matrix: None,
        }
    }

    /// Get id number
    pub fn get_id(&self) -> usize {
        let id_num: usize = self.id;
        id_num
    }

    /// Get the x-coords of nodes in Beam1D2N element
    pub fn xs(&self) -> [Dtype; 2] {
        let mut x_list = [0.0; 2];
        for i in 0..2 {
            x_list[i] = self.nodes[i].coord[0];
        }
        x_list
    }

    /// Get rod element length
    pub fn length(&self) -> Dtype {
        let x = self.xs();
        (x[0] - x[1]).abs()
    }

    /// Get nodes' displacement/rotate angle vector
    pub fn disps(&self) -> [Dtype; 4] {
        let mut disps = [0.0; 4];
        for idx in 0..2 {
            disps[2 * idx] = *self.nodes[idx].disps[0].borrow();
            disps[2 * idx + 1] = *self.nodes[idx].disps[1].borrow();
        }
        disps
    }

    /// Get nodes' force/moment vector
    pub fn forces(&self) -> [Dtype; 4] {
        let mut forces = [0.0; 4];
        for idx in 0..2 {
            forces[2 * idx] = *self.nodes[idx].forces[0].borrow();
            forces[2 * idx + 1] = *self.nodes[idx].forces[1].borrow();
        }
        forces
    }

    /// Shape matrix element N_i
    /// The shape mat of Beam1D2N element: [N1 N2 N3 N4]
    /// which is a 1x4 mat.    \xi = x/L
    /// N1 = 1*(1 + 0*\xi - 3*\xi^2 + 2*\xi^3)
    /// N2 = L*(0 + 1*\xi - 2*\xi^2 + 1*\xi^3)
    /// N3 = 1*(0 + 0*\xi + 3*\xi^2 - 2*\xi^3)
    /// N4 = L*(0 + 0*\xi - 1*\xi^2 + 1*\xi^3)
    fn shape_mat_i(&self, i: usize) -> impl Fn(Dtype) -> Dtype {
        // 参数i的取值范围：1～4
        let l = self.length();
        let l2 = l * l;
        let l3 = l * l * l;
        let a: [Dtype; 4] = [1.0, 0.0, 0.0, 0.0];
        let b: [Dtype; 4] = [0.0, l, 0.0, 0.0];
        let c: [Dtype; 4] = [-3.0, -2.0 * l, 3.0, -l];
        let d: [Dtype; 4] = [2.0, l, -2.0, l];
        move |x: Dtype| a[i] + b[i] * x / l + c[i] * x * x / l2 + d[i] * x * x * x / l3
    }

    /// Geometry matrix element B_i
    fn geometry_mat(&self, x: Dtype) -> SMatrix<Dtype, 1, 4> {
        let l = self.length();
        let l2 = l * l;
        let xi = x / l;
        SMatrix::<Dtype, 1, 4>::from([
            [(12.0 * xi - 6.0) / l2],
            [(6.0 * xi - 4.0) / l],
            [-(12.0 * xi - 6.0) / l2],
            [(6.0 * xi - 2.0) / l],
        ])
    }

    /// Calculate element stiffness matrix K
    /// arg moi: moment of inertia
    /// return a 4x4 matrix, elements are Dtype
    fn calc_k(&self, material_args: (Dtype, Dtype), moi: Dtype) -> [[Dtype; 4]; 4] {
        println!(
            "\n>>> Calculating Beam1D2N(#{})'s stiffness matrix k{} ......",
            self.id, self.id
        );

        let (ee, nu) = material_args;
        let l3 = self.length() * self.length() * self.length();
        let gauss_pt = (3.0 as Dtype).sqrt() / (6.0 as Dtype);
        let int_pts: [Dtype; 2] = [0.5 - gauss_pt, 0.5 + gauss_pt];

        let mut k_matrix = SMatrix::<Dtype, 4, 4>::from([[0.0; 4]; 4]);
        for idx in 0..2 {
            let b_mat = self.geometry_mat(int_pts[idx]);
            let core = b_mat.transpose() * b_mat;
            k_matrix += core;
        }

        let stiffness_matrix: [[Dtype; 4]; 4] = (ee * moi * k_matrix / l3).into();
        stiffness_matrix
    }
}
