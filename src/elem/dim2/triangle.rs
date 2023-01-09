use crate::{node::Node2D, Jacobian2D, K};
use na::*;
use std::fmt;

pub struct Tri2D3N<'tri> {
    pub id: usize,
    pub thick: f64,
    pub nodes: [&'tri Node2D; 3],
    pub k_matrix: Option<[[f64; 6]; 6]>,
}

impl<'tri> Tri2D3N<'tri> {
    /// Generate a 2D Tri2D3N element
    pub fn new(id: usize, thick: f64, nodes: [&Node2D; 3]) -> Tri2D3N {
        Tri2D3N {
            id,
            thick,
            nodes,
            k_matrix: None,
        }
    }

    /// Get the x-coords of nodes in tri element
    pub fn xs(&self) -> [f64; 3] {
        let mut x_list = [0.0; 3];
        for i in 0..3 {
            x_list[i] = self.nodes[i].coord[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in tri element
    pub fn ys(&self) -> [f64; 3] {
        let mut y_list = [0.0; 3];
        for i in 0..3 {
            y_list[i] = self.nodes[i].coord[1];
        }
        y_list
    }

    /// Get triangle element area value
    pub fn area(&self) -> f64 {
        let x = self.xs();
        let y = self.ys();
        let dx_21 = x[1] - x[0];
        let dx_31 = x[2] - x[0];
        let dy_21 = y[1] - y[0];
        let dy_31 = y[2] - y[0];
        0.5 * (dx_21 * dy_31 - dx_31 * dy_21).abs()
    }

    /// Calculate the Jacobian matrix of triangle element
    pub fn jacobian(&self) -> [[f64; 2]; 2] {
        let x: [f64; 3] = self.xs();
        let y: [f64; 3] = self.ys();
        let dx_21 = x[1] - x[0];
        let dx_31 = x[2] - x[0];
        let dy_21 = y[1] - y[0];
        let dy_31 = y[2] - y[0];
        [[dx_21, dx_31], [dy_21, dy_31]]
    }

    /// Get nodes' disps vector in tri element
    pub fn disps(&self) -> [f64; 6] {
        let mut disps = [0.0; 6];
        for idx in 0..3 {
            disps[2 * idx] = *self.nodes[idx].disps[0].borrow();
            disps[2 * idx + 1] = *self.nodes[idx].disps[1].borrow();
        }
        disps
    }

    /// Get nodes's force vector in tri element
    pub fn forces(&self) -> [f64; 6] {
        let mut forces = [0.0; 6];
        for idx in 0..3 {
            forces[2 * idx] = *self.nodes[idx].forces[0].borrow();
            forces[2 * idx + 1] = *self.nodes[idx].forces[1].borrow();
        }
        forces
    }

    /// Get any point's disps vector in tri element
    pub fn point_disp(&self, point_coord: [f64; 2]) -> [f64; 2] {
        let x = point_coord[0];
        let y = point_coord[1];
        let n0 = self.shape_mat_i(0usize)(x, y);
        let n1 = self.shape_mat_i(1usize)(x, y);
        let n2 = self.shape_mat_i(2usize)(x, y);

        let disps = self.disps();
        let u = n0 * disps[0] + n1 * disps[2] + n2 * disps[4];
        let v = n0 * disps[1] + n1 * disps[3] + n2 * disps[5];
        [u, v]
    }

    /// Get shape matrix element N_i
    fn shape_mat_i(&self, i: usize) -> impl Fn(f64, f64) -> f64 {
        /* The shape mat of tri elem:
         * |N1  0   N2  0   N3  0 |
         * |0   N1  0   N2  0   N3| */
        let area = self.area();
        let xs = self.xs();
        let ys = self.ys();
        let idx = |x: usize| (x % 3) as usize;
        let a = [
            xs[idx(0 + 1)] * ys[idx(0 + 2)] - xs[idx(0 + 2)] * ys[idx(0 + 1)],
            xs[idx(1 + 1)] * ys[idx(1 + 2)] - xs[idx(1 + 2)] * ys[idx(1 + 1)],
            xs[idx(2 + 1)] * ys[idx(2 + 2)] - xs[idx(2 + 2)] * ys[idx(2 + 1)],
        ];
        let b = [
            ys[idx(0 + 1)] - ys[idx(0 + 2)],
            ys[idx(1 + 1)] - ys[idx(1 + 2)],
            ys[idx(2 + 1)] - ys[idx(2 + 2)],
        ];
        let c = [
            xs[idx(0 + 2)] - xs[idx(0 + 1)],
            xs[idx(1 + 2)] - xs[idx(1 + 1)],
            xs[idx(2 + 2)] - xs[idx(2 + 1)],
        ];
        move |x: f64, y: f64| 0.5 * (a[i] + x * b[i] + y * c[i]) / area
    }

    /// Element's B matrix, B mat is the combination of diff(N)
    fn geometry_mat(&self, det_j: f64) -> SMatrix<f64, 3, 6> {
        let x: [f64; 3] = self.xs();
        let y: [f64; 3] = self.ys();
        let h_mat = SMatrix::<f64, 3, 4>::from([
            [y[2] - y[0], 0.0, x[0] - x[2]],
            [y[0] - y[1], 0.0, x[1] - x[0]],
            [0.0, x[0] - x[2], y[2] - y[0]],
            [0.0, x[1] - x[0], y[0] - y[1]],
        ]) / (det_j.abs());

        let q_mat = SMatrix::<f64, 4, 6>::from([
            [-1., -1., 0., 0.],
            [0., 0., -1., -1.],
            [1., 0., 0., 0.],
            [0., 0., 1., 0.],
            [0., 1., 0., 0.],
            [0., 0., 0., 1.],
        ]);

        h_mat * q_mat
    }

    /// Calculate element stiffness matrix K
    /// Return a 6x6 matrix, elements are f64
    fn calc_k(&self, material_args: (f64, f64)) -> [[f64; 6]; 6] {
        println!(
            "\n>>> Calculating Tri2D3N(#{})'s stiffness matrix k{} ......",
            self.id, self.id
        );
        let t = self.thick;
        let (ee, nu) = material_args;
        let elasticity_mat = SMatrix::<f64, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        let jacobian = Jacobian2D::from(self.jacobian());
        let det_j = jacobian.determinant();

        let b_mat = self.geometry_mat(det_j);
        // Gauss integration, area of standard tri is 0.5
        let core = b_mat.transpose() * elasticity_mat * b_mat * det_j;
        let stiffness_matrix: [[f64; 6]; 6] = (0.5 * t * core).into();
        stiffness_matrix
    }

    /// Get element's strain vector
    pub fn strain(&self) -> [f64; 3] {
        let jecobian_mat = Jacobian2D::from(self.jacobian());
        let b_mat = self.geometry_mat(jecobian_mat.determinant());
        let elem_nodes_disps = SMatrix::<f64, 6, 1>::from(self.disps());
        let strain_vlaue: [f64; 3] = (b_mat * elem_nodes_disps).into();
        strain_vlaue
    }

    /// Get element's strss vector
    pub fn stress(&self, material_args: (f64, f64)) -> [f64; 3] {
        let (ee, nu) = material_args;
        let elasticity_mat = SMatrix::<f64, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        let strain = SMatrix::<f64, 3, 1>::from(self.strain());
        let stress: [f64; 3] = (elasticity_mat * strain).into();
        stress
    }

    /// Print element's strain value
    pub fn print_strain(&self) {
        let strain = self.strain();
        println!(
            "\nelem[{}] strain:\n\tE_xx = {:-9.6}\n\tE_yy = {:-9.6}\n\tE_xy = {:-9.6}",
            self.id, strain[0], strain[1], strain[2]
        );
    }

    /// Print element's strss value
    pub fn print_stress(&self, material_args: (f64, f64)) {
        let stress = self.stress(material_args);
        println!(
            "\nelem[{}] stress:\n\tS_xx = {:-9.6}\n\tS_yy = {:-9.6}\n\tS_xy = {:-9.6}",
            self.id, stress[0], stress[1], stress[2]
        );
    }
}

impl<'tri> K for Tri2D3N<'tri> {
    type Kmatrix = [[f64; 6]; 6];

    /// Cache stiffness matrix for element
    fn k(&mut self, material: (f64, f64)) -> &Self::Kmatrix
    where
        Self::Kmatrix: std::ops::Index<usize>,
    {
        if self.k_matrix.is_none() {
            self.k_matrix.get_or_insert(self.calc_k(material))
        } else {
            self.k_matrix.as_ref().unwrap()
        }
    }

    /// Print element's stiffness matrix
    fn k_printer(&self, n_exp: f64) {
        if self.k_matrix.is_none() {
            println!(
                "!!! Tri2D3N#{}'s k mat is empty! call k() to calc it.",
                self.id
            );
        }

        print!("Tri2D3N k{} =  (* 10^{})\n[", self.id, n_exp as u8);
        for row in 0..6 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..6 {
                print!(
                    " {:>-10.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp))
                );
            }
            if row == 5 {
                println!("]]");
            } else {
                println!("]");
            }
        }
        println!("");
    }
}

impl fmt::Display for Tri2D3N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nElement_2D Info:\n\tId:     {}\n\tArea:   {}\n\tType:   Tri2D3N
\tNodes: {}\n\t       {}\n\t       {}",
            self.id,
            self.area(),
            self.nodes[0],
            self.nodes[1],
            self.nodes[2]
        )
    }
}
