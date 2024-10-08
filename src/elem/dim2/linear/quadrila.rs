use super::triangle::Tri2D3N;
use crate::{node::Node2D, Dtype, Jacobian2D, K};
use na::*;
use std::fmt::{self, Write};

pub struct Quad2D4N<'quad2d4n> {
    pub id: usize,
    pub thick: Dtype,
    pub nodes: [&'quad2d4n Node2D; 4],
    pub k_matrix: Option<[[Dtype; 8]; 8]>,
}

impl<'quad2d4n> Quad2D4N<'quad2d4n> {
    /// generate a new Rec2D4N element
    pub fn new(id: usize, thick: Dtype, nodes: [&Node2D; 4]) -> Quad2D4N {
        Quad2D4N {
            id,
            thick,
            nodes,
            k_matrix: None,
        }
    }

    /// Get id number
    pub fn get_id(&self) -> usize {
        let id_num: usize = self.id;
        id_num
    }

    /// get the x-coords of nodes in Rec2D4N element
    pub fn xs(&self) -> [Dtype; 4] {
        let mut x_list = [0.0; 4];
        for i in 0..4 {
            x_list[i] = self.nodes[i].coords[0];
        }
        x_list
    }

    /// get the y-coords of nodes in tri element
    pub fn ys(&self) -> [Dtype; 4] {
        let mut y_list = [0.0; 4];
        for i in 0..4 {
            y_list[i] = self.nodes[i].coords[1];
        }
        y_list
    }

    /// Get the rectangle element area
    pub fn area(&self) -> Dtype {
        let tri1: Tri2D3N = Tri2D3N {
            id: 0,
            thick: self.thick,
            nodes: [*self.nodes[0], *self.nodes[1], *self.nodes[2]],
            k_matrix: None,
        };
        let tri2: Tri2D3N = Tri2D3N {
            id: 0,
            thick: self.thick,
            nodes: [*self.nodes[3], *self.nodes[1], *self.nodes[2]],
            k_matrix: None,
        };
        tri1.area() + tri2.area()
    }

    /// Calculate the Jacobian matrix of quadrilateral
    pub fn jacobian(&self, xi_eta: [Dtype; 2]) -> [[Dtype; 2]; 2] {
        let x: [Dtype; 4] = self.xs();
        let y: [Dtype; 4] = self.ys();
        let xi: Dtype = xi_eta[0];
        let eta: Dtype = xi_eta[1];

        let p_xi: Dtype = 1.0 + xi;
        let n_xi: Dtype = 1.0 - xi;
        let p_eta: Dtype = 1.0 + eta;
        let n_eta: Dtype = 1.0 - eta;

        let x_d_xi: Dtype = 0.25 * (-n_eta * x[0] + n_eta * x[1] + p_eta * x[2] - p_eta * x[3]);
        let y_d_xi: Dtype = 0.25 * (-n_eta * y[0] + n_eta * y[1] + p_eta * y[2] - p_eta * y[3]);
        let x_d_eta: Dtype = 0.25 * (-n_xi * x[0] - p_xi * x[1] + p_xi * x[2] + n_xi * x[3]);
        let y_d_eta: Dtype = 0.25 * (-n_xi * y[0] - p_xi * y[1] + p_xi * y[2] + n_xi * y[3]);
        [[x_d_xi, y_d_xi], [x_d_eta, y_d_eta]]
    }

    /// Get nodes' disps vector in quad element
    pub fn disps(&self) -> [Dtype; 8] {
        let mut disps = [0.0; 8];
        for idx in 0..4 {
            disps[2 * idx] = self.nodes[idx].displs[0];
            disps[2 * idx + 1] = self.nodes[idx].displs[1];
        }
        disps
    }

    /// Get nodes' forces vector in quad element
    pub fn forces(&self) -> [Dtype; 8] {
        let mut forces = [0.0; 8];
        for idx in 0..4 {
            forces[2 * idx] = self.nodes[idx].forces[0];
            forces[2 * idx + 1] = self.nodes[idx].forces[1];
        }
        forces
    }

    /// Get any point's disps vector in quad element
    pub fn point_disp(&self, point_coord: [Dtype; 3]) -> [Dtype; 2] {
        let x = point_coord[0];
        let y = point_coord[1];
        let n0 = self.shape_mat_i(0usize)(x, y);
        let n1 = self.shape_mat_i(1usize)(x, y);
        let n2 = self.shape_mat_i(2usize)(x, y);
        let n3 = self.shape_mat_i(3usize)(x, y);

        let disps = self.disps();
        let u = n0 * disps[0] + n1 * disps[2] + n2 * disps[4] + n3 * disps[6];
        let v = n0 * disps[1] + n1 * disps[3] + n2 * disps[5] + n3 * disps[7];
        [u, v]
    }

    /// Get shape matrix element N_i
    fn shape_mat_i(&self, i: usize) -> impl Fn(Dtype, Dtype) -> Dtype {
        /* The shape mat of quad elem:
         * |N1   0    N2   0    N3   0    N4   0 |
         * |0    N1   0    N2   0    N3   0    N4|  */
        let x_sign = [-1.0, 1.0, 1.0, -1.0];
        let y_sign = [-1.0, -1.0, 1.0, 1.0];
        move |x: Dtype, y: Dtype| {
            0.25 * (1.0 + x_sign[i] * x + y_sign[i] * y + x_sign[i] * y_sign[i] * x * y)
        }
    }

    /// Strain in x-y coord,
    /// epsilon = Du/Dx = Du/D(xi) * D(xi)/Dx = Du/D(xi) * J^(-1)
    /// return conmbination of D(xi)/Dx which is h_mat
    fn h_mat(&self, jacobian: [[Dtype; 2]; 2], det_j: Dtype) -> SMatrix<Dtype, 3, 4> {
        SMatrix::<Dtype, 3, 4>::from([
            [jacobian[1][1], 0.0, -jacobian[1][0]],
            [-jacobian[0][1], 0.0, jacobian[0][0]],
            [0.0, -jacobian[1][0], jacobian[1][1]],
            [0.0, jacobian[0][0], -jacobian[0][1]],
        ]) / det_j
    }

    /// Strain in xi-eta coord,
    /// epsilon = Du/D(xi) = D(N*q)/D(xi) = [D(N)/D(xi)]*q
    /// return [D(N)/D(xi)] which is q_mat = diff(N(xi, eta))
    fn q_mat(&self, xi_eta: [Dtype; 2]) -> SMatrix<Dtype, 4, 8> {
        let p_xi: Dtype = 1.0 + xi_eta[0];
        let n_xi: Dtype = 1.0 - xi_eta[0];
        let p_eta: Dtype = 1.0 + xi_eta[1];
        let n_eta: Dtype = 1.0 - xi_eta[1];

        SMatrix::<Dtype, 4, 8>::from([
            [-n_eta, -n_xi, 0.0, 0.0],
            [0.0, 0.0, -n_eta, -n_xi],
            [n_eta, -p_xi, 0.0, 0.0],
            [0.0, 0.0, n_eta, -p_xi],
            [p_eta, p_xi, 0.0, 0.0],
            [0.0, 0.0, p_eta, p_xi],
            [-p_eta, n_xi, 0.0, 0.0],
            [0.0, 0.0, -p_eta, n_xi],
        ]) * 0.25
    }

    /// Element's B matrix which is the conmbination of diff(N(x,y))
    fn geometry_mat(
        &self,
        j_raw: [[Dtype; 2]; 2],
        det_j: Dtype,
        xi_eta: [Dtype; 2],
    ) -> SMatrix<Dtype, 3, 8> {
        self.h_mat(j_raw, det_j) * self.q_mat(xi_eta)
    }

    /// Calculate element stiffness matrix K
    /// return a 8x8 matrix, elements are Dtype
    fn calc_k(&self, material_args: (Dtype, Dtype)) -> [[Dtype; 8]; 8] {
        println!(
            "\n>>> Calculating Quad2D4N(#{})'s stiffness matrix k{} ......",
            self.id, self.id
        );
        let (ee, nu) = material_args;
        let elasticity_mat = SMatrix::<Dtype, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        // 4 point Gauss integration, area of standard rectangle is 4.0
        let gauss_pt = (1.0 as Dtype) / ((3.0 as Dtype).sqrt());
        let int_pts = [
            [[-gauss_pt, -gauss_pt], [gauss_pt, -gauss_pt]],
            [[gauss_pt, gauss_pt], [-gauss_pt, gauss_pt]],
        ];

        let mut k_matrix = SMatrix::<Dtype, 8, 8>::from([[0.0; 8]; 8]);
        // 2x2的4点Gauss积分，4个点权重都是1
        for row in 0..2 {
            for col in 0..2 {
                let j_raw = self.jacobian(int_pts[row][col]);
                let j = Jacobian2D::from(j_raw);
                let det_j = j.determinant().abs();

                let b_mat = self.geometry_mat(j_raw, det_j, int_pts[row][col]);
                let core = b_mat.transpose() * elasticity_mat * b_mat * det_j;
                k_matrix += core;
            }
        }

        let stiffness_matrix: [[Dtype; 8]; 8] = (self.thick * k_matrix).into();
        stiffness_matrix
    }

    /// Get element's strain vector
    fn calc_strain(&self, xi_eta: [Dtype; 2]) -> [Dtype; 3] {
        let j_raw = self.jacobian(xi_eta);
        let det_j = Jacobian2D::from(self.jacobian(xi_eta)).determinant();
        let b_mat = self.geometry_mat(j_raw, det_j, xi_eta);
        let elem_nodes_disps = SMatrix::<Dtype, 8, 1>::from(self.disps());
        let strain_vector: [Dtype; 3] = (b_mat * elem_nodes_disps).into();
        strain_vector
    }

    /// Get element's strss vector
    fn calc_stress(&self, xi_eta: [Dtype; 2], material_args: (Dtype, Dtype)) -> [Dtype; 3] {
        let (ee, nu) = material_args;
        let elasticity_mat = SMatrix::<Dtype, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        let strain = SMatrix::<Dtype, 3, 1>::from(self.calc_strain(xi_eta));
        let stress: [Dtype; 3] = (elasticity_mat * strain).into();
        stress
    }

    /// Print element's strain value
    pub fn print_strain(&self, xi_eta: [Dtype; 2]) {
        let strain = self.calc_strain(xi_eta);
        println!(
            "\nelem[{}] strain:\n\tE_xx = {:-9.6}\n\tE_yy = {:-9.6}\n\tE_xy = {:-9.6}",
            self.id, strain[0], strain[1], strain[2]
        );
    }

    /// Print element's stress value
    pub fn print_stress(&self, xi_eta: [Dtype; 2], material_args: (Dtype, Dtype)) {
        let stress = self.calc_stress(xi_eta, material_args);
        println!(
            "\nelem[{}] stress:\n\tS_xx = {:-9.6}\n\tS_yy = {:-9.6}\n\tS_xy = {:-9.6}",
            self.id, stress[0], stress[1], stress[2]
        );
    }
}

/// Implement zhm::K trait for quad element
impl<'quad2d4n> K for Quad2D4N<'quad2d4n> {
    type Kmatrix = [[Dtype; 8]; 8];

    /// Cache stiffness matrix for quad element
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

    /// Print quad element's stiffness matrix
    fn k_printer(&self, n_exp: Dtype) {
        if self.k_matrix.is_none() {
            panic!(
                "!!! Quad2D4N#{}'s k mat is empty! call k() to calc it.",
                self.id
            )
        }

        print!("Quad2D4N k{} =  (* 10^{})\n[", self.id, n_exp as i32);
        for row in 0..8 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..8 {
                print!(
                    " {:>-10.6}",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                );
            }
            if row == 7 {
                println!("]]");
            } else {
                println!("]");
            }
        }
        print!("\n");
    }

    /// Return quad elem's stiffness matrix's format string
    fn k_string(&self, n_exp: Dtype) -> String {
        if self.k_matrix.is_none() {
            panic!(
                "!!! Quad2D4N#{}'s k mat is empty! call k() to calc it.",
                self.id
            )
        }

        let mut k_matrix = String::new();
        for row in 0..8 {
            if row == 0 {
                write!(k_matrix, "[[").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, " [").expect("!!! Write tri k_mat failed!");
            }
            for col in 0..8 {
                write!(
                    k_matrix,
                    " {:>-10.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                )
                .expect("!!! Write tri k_mat failed!");
            }
            if row == 5 {
                write!(k_matrix, "]]").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, "]\n").expect("!!! Write tri k_mat failed!");
            }
        }
        k_matrix
    }

    /// Get the strain at (xi, eta) inside the element
    fn strain(&self, xyz: [Dtype; 3]) -> Vec<Dtype> {
        self.calc_strain([xyz[0], xyz[1]]).to_vec()
    }

    /// Get the stress at (xi, eta) inside the element
    fn stress(&self, xyz: [Dtype; 3], material: (Dtype, Dtype)) -> Vec<Dtype> {
        self.calc_stress([xyz[0], xyz[1]], material).to_vec()
    }

    /// Get element's info string
    fn info(&self) -> String {
        format!("\n--------------------------------------------------------------------\nElement_2D Info:\n\tId:     {}\n\tArea:   {}\n\tType:   Quad2D4N
\tNodes: {}\n\t       {}\n\t       {}\n\t       {}\n\t",
            self.id,
            self.area(),
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
            self.nodes[3],
        )
    }

    /// Get element's id number
    fn id(&self) -> usize {
        self.id
    }
}

impl fmt::Display for Quad2D4N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nElement_2D Info:\n\tId:     {}\n\tArea:   {}\n\tType:   Rec2D4N
\tNodes: {}\n\t       {}\n\t       {}\n\t       {}",
            self.id,
            self.area(),
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
            self.nodes[3]
        )
    }
}
