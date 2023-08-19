use crate::{node::Node2D, Dtype, K};
use na::*;
use std::fmt::{self, Write};

/// 细长梁以及纯弯梁，抽象成具有如下特征的理论模型：
///   1) 仅用x坐标刻画所有位置(即细长之含义)；
///   2) 变形主要为垂直于x轴的挠度
///   3) 节点受载荷为横力和弯矩
/// 但是，节点的自由度不为1，包含：挠度(y)、转角(theta)，
/// 此亦是梁单元为非协调单元的原因。
pub struct Beam1D2N<'beam1d2n> {
    pub id: usize,
    pub moi: Dtype, // moment of inertia
    pub sec_area: Dtype,
    pub nodes: [&'beam1d2n Node2D; 2],
    pub k_matrix: Option<[[Dtype; 4]; 4]>,
}

impl<'beam1d2n> Beam1D2N<'beam1d2n> {
    /// Generate a Beam1D2N element
    pub fn new(id: usize, moi: Dtype, sec_area: Dtype, nodes: [&Node2D; 2]) -> Beam1D2N {
        Beam1D2N {
            id,
            moi,
            sec_area,
            nodes,
            k_matrix: None,
        }
    }

    /// Get id number
    pub fn get_id(&self) -> usize {
        let id_num: usize = self.id.clone();
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

    /// Get point(in beam element)'s disp
    /// 用节点的位移插值得到单元内任意点的位移，由于是纯弯梁
    /// 这里的位移指每个节点的纵向位移v和转角\theta
    pub fn point_disp(&self, point_coord: [Dtype; 1]) -> [Dtype; 2] {
        let x = point_coord[0];
        let n0 = self.shape_mat_i(0usize)(x);
        let n1 = self.shape_mat_i(1usize)(x);
        let n2 = self.shape_mat_i(2usize)(x);
        let n3 = self.shape_mat_i(3usize)(x);
        let disps = self.disps();
        let v = n0 * disps[0] + n2 * disps[2];
        let theta = n1 * disps[1] + n3 * disps[3];
        [v, theta]
    }

    /// Shape matrix element N_i
    /// The shape mat of Beam1D2N element: [N1 N2 N3 N4]
    /// which is a 1x4 mat.    \xi = x/L
    /// N1 = 1*(1 + 0*\xi - 3*\xi^2 + 2*\xi^3)
    /// N2 = L*(0 + 1*\xi - 2*\xi^2 + 1*\xi^3)
    /// N3 = 1*(0 + 0*\xi + 3*\xi^2 - 2*\xi^3)
    /// N4 = L*(0 + 0*\xi - 1*\xi^2 + 1*\xi^3)
    fn shape_mat_i(&self, i: usize) -> impl Fn(Dtype) -> Dtype {
        // 参数i的取值范围：0～3
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
    fn geometry_mat(&self, xi: Dtype) -> SMatrix<Dtype, 1, 4> {
        let l = self.length();
        SMatrix::<Dtype, 1, 4>::from([
            [12.0 * xi - 6.0],
            [(6.0 * xi - 4.0) * l],
            [-12.0 * xi + 6.0],
            [(6.0 * xi - 2.0) * l],
        ])
    }

    /// Calculate element stiffness matrix K
    /// arg moi: moment of inertia
    /// return a 4x4 matrix, elements are Dtype
    pub fn calc_k(&self, material_args: (Dtype, Dtype)) -> [[Dtype; 4]; 4] {
        println!(
            "\n>>> Calculating Beam1D2N(#{})'s stiffness matrix k{} ......",
            self.id, self.id
        );

        let moi = self.moi;
        let (ee, _nu) = material_args;
        let l = self.length();
        let l3 = l * l * l;
        let gauss_pt: Dtype = (3.0).sqrt() / 6.0;
        let int_pts: [Dtype; 2] = [0.5 - gauss_pt, 0.5 + gauss_pt];

        let mut k_matrix = SMatrix::<Dtype, 4, 4>::from([[0.0; 4]; 4]);
        for idx in 0..2 {
            let b_mat = self.geometry_mat(int_pts[idx]);
            let core = 0.5 * b_mat.transpose() * b_mat;
            k_matrix = k_matrix + core;
        }

        let stiffness_matrix: [[Dtype; 4]; 4] = (ee * moi * k_matrix / l3).into();
        stiffness_matrix
    }

    /// Get the strain at point (x) inside the element
    pub fn calc_strain(&self, xy: [Dtype; 2]) -> [Dtype; 1] {
        let y = xy[1];
        let xi = xy[0] / self.length();
        let b_mat = self.geometry_mat(xi);
        let elem_node_disps = SMatrix::<Dtype, 4, 1>::from(self.disps());
        let strain_vector: [Dtype; 1] = (-y * b_mat * elem_node_disps).into();
        strain_vector
    }

    /// Get the stress at point (x) inside the element
    pub fn calc_stress(&self, xy: [Dtype; 2], material_args: (Dtype, Dtype)) -> [Dtype; 1] {
        let (ee, _nu) = material_args;
        [ee * self.calc_strain(xy)[0]]
    }
}

/// Implement zhm::K trait for beam1d2n element
impl<'beam1d2n> K for Beam1D2N<'beam1d2n> {
    type Kmatrix = [[Dtype; 4]; 4];

    /// Cache stiffness matrix for beam1d2n element
    fn k(&mut self, material: (Dtype, Dtype)) -> &Self::Kmatrix
    where
        Self::Kmatrix: std::ops::Index<usize>,
    {
        if self.k_matrix.is_none() {
            self.k_matrix.get_or_insert(self.calc_k(material))
        } else {
            self.k_matrix.as_ref().unwrap()
        }
    }

    /// Print beam1d2n element's stiffness matrix
    fn k_printer(&self, n_exp: Dtype) {
        if self.k_matrix.is_none() {
            panic!(
                "!!! Beam1D2N#{}'s k mat is empty! call k() to calc it.",
                self.id
            );
        }

        print!("\n Beam1D2N k{} = (* 10^{})\n[", self.id, n_exp as u8);
        for row in 0..4 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..4 {
                print!(
                    " {:>-10.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                );
            }
            if row == 3 {
                print!("]]");
            } else {
                print!("]");
            }
        }
        print!("\n");
    }

    /// Return Beam1D2N elem's stiffness matrix's format string
    fn k_string(&self, n_exp: Dtype) -> String {
        let mut k_matrix = String::new();
        for row in 0..4 {
            if row == 0 {
                write!(k_matrix, "[[").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, " [").expect("!!! Write tri k_mat failed!");
            }
            for col in 0..4 {
                write!(
                    k_matrix,
                    " {:>-10.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                )
                .expect("!!! Write element's stiffness mat string failed.");
            }
            if row == 3 {
                write!(k_matrix, "]]").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, "]\n").expect("!!! Write tri k_mat failed!");
            }
        }
        k_matrix
    }

    /// Get the strain at point (x) inside the element
    fn strain(&self, xy: [Dtype; 3]) -> Vec<Dtype> {
        self.calc_strain([xy[0], xy[1]]).to_vec()
    }

    /// Get the stress at point (x) inside the element
    fn stress(&self, xy: [Dtype; 3], material_args: (Dtype, Dtype)) -> Vec<Dtype> {
        self.calc_stress([xy[0], xy[1]], material_args).to_vec()
    }

    /// Get element's info string
    fn info(&self) -> String {
        format!(
            "\n--------------------------------------------------------------------\nElement_1D Info:\n\tId:     {}\n\tCrossArea:   {}\n\tType:   Beam1D2N
\tNodes: {}\n\t       {}\n\t",self.id, self.sec_area, self.nodes[0], self.nodes[1],
        )
    }

    /// Get element's id number
    fn id(&self) -> usize {
        self.id
    }
}

impl fmt::Display for Beam1D2N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nElement_1D Info:\n\tId:     {}\n\tArea:   {}\n\tType:   Beam1D2N (Euler beam) 
\tNodes: {}\n\t       {}\n\t",
            self.id, self.sec_area, self.nodes[0], self.nodes[1],
        )
    }
}
