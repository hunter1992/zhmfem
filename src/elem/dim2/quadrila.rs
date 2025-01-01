use super::triangle::Tri2D3N;
use crate::{node::Node2D, Dtype, Jacobian2D, K};
use na::SMatrix;
use std::fmt::{self, Write};

pub struct Quad2D4N<'quad2d4n> {
    pub id: usize,
    pub thick: Dtype,
    pub nodes: [&'quad2d4n Node2D; 4],
    pub k_matrix: Option<[[Dtype; 8]; 8]>,
    pub material: &'quad2d4n (Dtype, Dtype),
}

impl<'quad2d4n> Quad2D4N<'quad2d4n> {
    /// Generate a new Quad2D4N element
    pub fn new(
        id: usize,
        thick: Dtype,
        nodes: [&'quad2d4n Node2D; 4],
        material: &'quad2d4n (Dtype, Dtype),
    ) -> Self {
        Quad2D4N {
            id,
            thick,
            nodes,
            k_matrix: None,
            material,
        }
    }

    /// Set element material_args
    pub fn set_material(&mut self, material_args: &'quad2d4n (Dtype, Dtype)) {
        self.material = material_args;
    }

    /// Get the rectangle element area
    pub fn area(&self) -> Dtype {
        let tri1: Tri2D3N = Tri2D3N {
            id: 0,
            thick: self.thick,
            nodes: [self.nodes[0], self.nodes[1], self.nodes[2]],
            k_matrix: None,
            material: self.material,
        };
        let tri2: Tri2D3N = Tri2D3N {
            id: 1,
            thick: self.thick,
            nodes: [self.nodes[2], self.nodes[3], self.nodes[0]],
            k_matrix: None,
            material: self.material,
        };
        tri1.area() + tri2.area()
    }

    /// Get the x-coords of nodes in tri element
    pub fn get_nodes_xcoords(&self) -> [Dtype; 4] {
        let mut x_list = [0.0; 4];
        for i in 0..4 {
            x_list[i] = self.nodes[i].coords[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in tri element
    pub fn get_nodes_ycoords(&self) -> [Dtype; 4] {
        let mut y_list = [0.0; 4];
        for i in 0..4 {
            y_list[i] = self.nodes[i].coords[1];
        }
        y_list
    }

    /// Get nodes' disps vector in tri element
    pub fn get_nodes_displacement(&self) -> [Dtype; 8] {
        let mut disps = [0.0; 8];
        for idx in 0..4 {
            disps[2 * idx] = self.nodes[idx].displs.borrow()[0];
            disps[2 * idx + 1] = self.nodes[idx].displs.borrow()[1];
        }
        disps
    }

    /// Get nodes's force vector in tri element
    pub fn get_nodes_force(&self) -> [Dtype; 8] {
        let mut forces = [0.0; 8];
        for idx in 0..4 {
            forces[2 * idx] = self.nodes[idx].forces.borrow()[0];
            forces[2 * idx + 1] = self.nodes[idx].forces.borrow()[1];
        }
        forces
    }

    /// Get shape matrix element N_i
    ///
    /// -------------------
    ///   i  |  xi  |  eta
    /// -------------------
    ///   1  |  -1  |  -1   (从第三象限开始)
    ///   2  |   1  |  -1
    ///   3  |   1  |   1
    ///   4  |  -1  |   1
    /// -------------------
    ///
    /// N_i(xi, eta) = 0.25 * (1 + xi_i * xi) * (1 + eta_i * eta)
    fn shape_mat_i(&self, ith: usize) -> impl Fn(Dtype, Dtype) -> Dtype {
        /* The shape mat of quad elem:
         * |N1   0    N2   0    N3   0    N4   0 |
         * |0    N1   0    N2   0    N3   0    N4|  */
        let x_sign = [-1.0, 1.0, 1.0, -1.0];
        let y_sign = [-1.0, -1.0, 1.0, 1.0];
        move |x: Dtype, y: Dtype| {
            0.25 * (1.0 + x_sign[ith] * x + y_sign[ith] * y + x_sign[ith] * y_sign[ith] * x * y)
        }
    }

    /// Get any point's disps vector in quad element
    pub fn point_disp(&self, xi_eta: [Dtype; 3]) -> [Dtype; 2] {
        let s = xi_eta[0];
        let t = xi_eta[1];
        let n0 = self.shape_mat_i(0usize)(s, t);
        println!(">>> Called <<<");
        let n1 = self.shape_mat_i(1usize)(s, t);
        let n2 = self.shape_mat_i(2usize)(s, t);
        let n3 = self.shape_mat_i(3usize)(s, t);

        let disps = self.get_nodes_displacement();
        let u = n0 * disps[0] + n1 * disps[2] + n2 * disps[4] + n3 * disps[6];
        let v = n0 * disps[1] + n1 * disps[3] + n2 * disps[5] + n3 * disps[7];
        [u, v]
    }

    /// Calculate the Jacobian matrix of quadrilateral
    fn jacobian(&self, xi_eta: [Dtype; 3]) -> [[Dtype; 2]; 2] {
        let x: [Dtype; 4] = self.get_nodes_xcoords();
        let y: [Dtype; 4] = self.get_nodes_ycoords();
        let s: Dtype = xi_eta[0];
        let t: Dtype = xi_eta[1];

        let p_s: Dtype = 1.0 + s; // positive xi
        let n_s: Dtype = 1.0 - s;
        let p_t: Dtype = 1.0 + t;
        let n_t: Dtype = 1.0 - t;

        let dx_ds: Dtype = 0.25 * (-n_t * x[0] + n_t * x[1] + p_t * x[2] - p_t * x[3]);
        let dy_ds: Dtype = 0.25 * (-n_t * y[0] + n_t * y[1] + p_t * y[2] - p_t * y[3]);
        let dx_dt: Dtype = 0.25 * (-n_s * x[0] - p_s * x[1] + p_s * x[2] + n_s * x[3]);
        let dy_dt: Dtype = 0.25 * (-n_s * y[0] - p_s * y[1] + p_s * y[2] + n_s * y[3]);
        [[dx_ds, dy_ds], [dx_dt, dy_dt]]
    }

    /// Calculate H matrix
    /// The H matrix transforms the strain in the parameter coordinate system to
    /// the physical coordinate system
    ///
    /// epsilon = Du/Dx
    ///         = Du/Ds * Ds/Dx
    ///         = J^(-1) * Du/Ds
    ///         = H * Epsilon
    ///
    /// |     du/dx     |        | du/ds |
    /// |     dv/dy     |  = H * | du/dt |
    /// | du/dy + dv/dx |        | dv/ds |
    ///                          | dv/dt |
    ///
    /// epsilon: Strain in physical coordinate system
    /// Epsilon: Strain in parametric coordinate system
    fn geometry_transform_mat(&self, xi_eta: [Dtype; 3]) -> SMatrix<Dtype, 3, 4> {
        let j_mat = self.jacobian(xi_eta);
        let det_j: Dtype = j_mat[0][0] * j_mat[1][1] - j_mat[0][1] * j_mat[1][0];
        if det_j < 0.0 {
            panic!("!!! Det of Jacobian < 0 !");
        }
        let h = SMatrix::<Dtype, 3, 4>::from([
            [j_mat[1][1], 0.0, -j_mat[1][0]],
            [-j_mat[0][1], 0.0, j_mat[0][0]],
            [0.0, -j_mat[1][0], j_mat[1][1]],
            [0.0, j_mat[0][0], -j_mat[0][1]],
        ]);
        h / det_j
    }

    /// Calculate Q matrix
    /// Q matrix is the  differentiation of N(s, t):
    ///
    /// | du/ds |   | d /ds,   0   |   | u |
    /// | du/dt | = | d /dt,   0   | * |   |
    /// | dv/ds |   |   0  , d /ds |   | v |
    /// | dv/dt |   |   0  , d /dt |   |   |
    ///                                                                     | u1 |
    ///             | d /ds,   0   |                                        | v1 |
    ///           = | d /dt,   0   | * | N1, 0 , N2, 0 , N3, 0 , N4, 0  | * | u2 |
    ///             |   0  , d /ds |   | 0 , N1, 0 , N2, 0 , N3, 0 , N4 |   | v2 |
    ///             |   0  , d /dt |                                        | u3 |
    ///                                                                     | v3 |
    ///           = Q * [ u1 v1 u2 v2 u3 v3 u4 v4 ]^T                       | u4 |
    ///
    /// Q = | d /ds,   0   |
    ///     | d /dt,   0   | * | N1, 0 , N2, 0 , N3, 0 , N4, 0  |
    ///     |   0  , d /ds |   | 0 , N1, 0 , N2, 0 , N3, 0 , N4 |
    ///     |   0  , d /dt |
    ///
    ///   = | (1 + eta_1 * eta) * xi_1 ... |
    ///     | .                            |       *  0.25
    ///     | .                            |
    ///     | ...                          |(4x8)
    fn diff_shape_mat(&self, xi_eta: [Dtype; 3]) -> SMatrix<Dtype, 4, 8> {
        let s: Dtype = xi_eta[0];
        let t: Dtype = xi_eta[1];
        let p_s: Dtype = 1.0 + s;
        let n_s: Dtype = 1.0 - s;
        let p_t: Dtype = 1.0 + t;
        let n_t: Dtype = 1.0 - t;

        SMatrix::<Dtype, 4, 8>::from([
            [-n_t, -n_s, 0.0, 0.0],
            [0.0, 0.0, -n_t, -n_s],
            [n_t, -p_s, 0.0, 0.0],
            [0.0, 0.0, n_t, -p_s],
            [p_t, p_s, 0.0, 0.0],
            [0.0, 0.0, p_t, p_s],
            [-p_t, n_s, 0.0, 0.0],
            [0.0, 0.0, -p_t, n_s],
        ]) * 0.25
    }

    /// Element's B matrix which is the conmbination of diff(N(x,y))
    fn geometry_mat(&self, xi_eta: [Dtype; 3]) -> SMatrix<Dtype, 3, 8> {
        let h_mat = self.geometry_transform_mat(xi_eta);
        let q_mat = self.diff_shape_mat(xi_eta);
        h_mat * q_mat
    }

    /// Calculate element stiffness matrix K
    /// return a 8x8 matrix, elements are Dtype
    fn calc_k(&self) -> [[Dtype; 8]; 8] {
        println!(
            "\n>>> Calculating Quad2D4N(#{})'s stiffness matrix k{} ......",
            self.id, self.id
        );
        let (ee, nu) = *self.material;
        let elasticity_mat = SMatrix::<Dtype, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        // 4 point Gauss integration, area of standard rectangle is 4.0
        let gauss_pt = (1.0 as Dtype) / ((3.0 as Dtype).sqrt());
        let int_pts = [
            [[-gauss_pt, -gauss_pt, 0.0], [gauss_pt, -gauss_pt, 0.0]],
            [[gauss_pt, gauss_pt, 0.0], [-gauss_pt, gauss_pt, 0.0]],
        ];

        let mut k_matrix = SMatrix::<Dtype, 8, 8>::from([[0.0; 8]; 8]);
        // 2x2的4点Gauss积分，4个点权重都是1
        for row in 0..2 {
            for col in 0..2 {
                let j_raw = self.jacobian(int_pts[row][col]);
                let j = Jacobian2D::from(j_raw);
                let det_j = j.determinant().abs();

                let b_mat = self.geometry_mat(int_pts[row][col]);
                let core = b_mat.transpose() * elasticity_mat * b_mat * det_j;
                k_matrix += core;
            }
        }

        let stiffness_matrix: [[Dtype; 8]; 8] = (self.thick * k_matrix).into();
        stiffness_matrix
    }

    /// Get element's strain vector
    fn calc_strain(&self, xi_eta: [Dtype; 3]) -> [Dtype; 3] {
        let b_mat = self.geometry_mat(xi_eta);
        let elem_nodes_disps = SMatrix::<Dtype, 8, 1>::from(self.get_nodes_displacement());
        let strain_vector: [Dtype; 3] = (b_mat * elem_nodes_disps).into();
        strain_vector
    }

    /// Get element's strss vector
    fn calc_stress(&self, xi_eta: [Dtype; 3]) -> [Dtype; 3] {
        let (ee, nu) = *self.material;
        let elasticity_mat = SMatrix::<Dtype, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        let strain = SMatrix::<Dtype, 3, 1>::from(self.calc_strain(xi_eta));
        let stress: [Dtype; 3] = (elasticity_mat * strain).into();
        stress
    }

    /// Get element node's stress
    pub fn nodes_stress(&self) -> [[Dtype; 3]; 4] {
        let mut sigma: [[Dtype; 3]; 4] = [[0.0; 3]; 4];
        let nodes_coords: [[Dtype; 3]; 4] = [
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
            [-1.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
        ];
        for idx in 0..4 {
            sigma[idx] = self.calc_stress(nodes_coords[idx]);
        }
        sigma
    }

    /// Print element's strain value
    pub fn print_strain(&self, xi_eta: [Dtype; 3]) {
        let strain = self.calc_strain(xi_eta);
        println!(
            "\nelem[{}] strain:\n\tE_xx = {:-16.6}\n\tE_yy = {:-16.6}\n\tE_xy = {:-16.6}",
            self.id, strain[0], strain[1], strain[2]
        );
    }

    /// Print element's stress value
    pub fn print_stress(&self, xi_eta: [Dtype; 3]) {
        let stress = self.calc_stress(xi_eta);
        println!(
            "\nelem[{}] stress:\n\tS_xx = {:-16.6}\n\tS_yy = {:-16.6}\n\tS_xy = {:-16.6}",
            self.id, stress[0], stress[1], stress[2]
        );
    }
}

/// Implement zhm::K trait for quadrilateral element
impl<'quad2d4n> K for Quad2D4N<'quad2d4n> {
    type Kmatrix = [[Dtype; 8]; 8];

    /// Cache stiffness matrix for quad element
    fn k(&mut self) -> &Self::Kmatrix
    where
        Self::Kmatrix: std::ops::Index<usize>,
    {
        if self.k_matrix.is_none() {
            self.k_matrix.get_or_insert(self.calc_k())
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
                    " {:>-12.6}",
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
                    " {:>-12.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                )
                .expect("!!! Write tri k_mat failed!");
            }
            if row == 7 {
                write!(k_matrix, "]]").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, "]\n").expect("!!! Write tri k_mat failed!");
            }
        }
        k_matrix
    }

    /// Get the strain at (xi, eta) inside the element
    fn strain(&self, xyz: [Dtype; 3]) -> Vec<Dtype> {
        self.calc_strain(xyz).to_vec()
    }

    /// Get the stress at (xi, eta) inside the element
    fn stress(&self, xyz: [Dtype; 3]) -> Vec<Dtype> {
        self.calc_stress(xyz).to_vec()
    }

    /// Get element's info string
    fn info(&self, n_exp: Dtype) -> String {
        let gauss_pt = (1.0 as Dtype) / ((3.0 as Dtype).sqrt());
        let int_pts: [[[Dtype; 3]; 2]; 2] = [
            [[-gauss_pt, -gauss_pt, 0.0], [gauss_pt, -gauss_pt, 0.0]],
            [[gauss_pt, gauss_pt, 0.0], [-gauss_pt, gauss_pt, 0.0]],
        ];
        format!("\n-----------------------------------------------------------------------------\nElem_Quad2D4N:\n\tId:\t{}\n\tArea: {:-12.6}\n\tMats: {:-12.6} (Young's modulus)\n\t      {:-12.6} (Poisson's ratio)\n\tNodes:{}{}{}{}\n\tStrain:\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\n\tStress:\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\n\tStiffness Matrix K{} =  (*10^{})\n{}",
            self.id,
            self.area(),
            self.material.0,
            self.material.1,
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
            self.nodes[3],
            self.strain(int_pts[0][0]),
            self.strain(int_pts[0][1]),
            self.strain(int_pts[1][0]),
            self.strain(int_pts[1][1]),
            self.stress(int_pts[0][0]),
            self.stress(int_pts[0][1]),
            self.stress(int_pts[1][0]),
            self.stress(int_pts[1][1]),
            self.id(),
            n_exp,
            self.k_string(n_exp),
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
            "\nElem_Quad2D4N:\n\tId:\t{}\n\tArea: {:-12.6}\n\tMats: {:-12.6} (Young's modulus)\n\t      {:-12.6} (Poisson's ratio)\n\tNodes:{}{}{}{}",
            self.id,
            self.area(),
            self.material.0,
            self.material.1,
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
            self.nodes[3],
        )
    }
}
