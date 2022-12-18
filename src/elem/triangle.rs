extern crate nalgebra as na;

use crate::{node::Node2D, K};
use na::*;
use std::fmt;

type Jacobian2x2f = SMatrix<f64, 2, 2>;

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

    /// Get the x coords of nodes in tri element
    pub fn xs(&self) -> [f64; 3] {
        let mut x_list = [0.0; 3];
        for i in 0..3 {
            x_list[i] = self.nodes[i].coord[0];
        }
        x_list
    }

    /// Get the y coords of nodes in tri element
    pub fn ys(&self) -> [f64; 3] {
        let mut y_list = [0.0; 3];
        for i in 0..3 {
            y_list[i] = self.nodes[i].coord[1];
        }
        y_list
    }

    /// Get nodes's disps in a element
    pub fn disps(&self) -> [f64; 6] {
        let mut disps = [0.0; 6];
        for idx in 0..3 {
            disps[2 * idx] = *self.nodes[idx].disps[0].borrow();
            disps[2 * idx + 1] = *self.nodes[idx].disps[1].borrow();
        }
        disps
    }

    /// Calculate the Jacobian matrix of triangle element
    pub fn jacobian(&self) -> [[f64; 2]; 2] {
        let x: [f64; 3] = self.xs();
        let y: [f64; 3] = self.ys();
        let dx21 = x[1] - x[0];
        let dx31 = x[2] - x[0];
        let dy21 = y[1] - y[0];
        let dy31 = y[2] - y[0];
        [[dx21, dx31], [dy21, dy31]]
    }

    /// Calculate element stiffness matrix K
    /// Return a 6x6 matrix, elements are f64
    fn calc_k(&self, material_args: (f64, f64)) -> [[f64; 6]; 6] {
        println!(
            "\n>>> Calculating Tri2D3N#{}'s stiffness matrix k{} ......",
            self.id, self.id
        );
        let t = self.thick;
        let (ee, nu) = material_args;
        let elasticity_mat = SMatrix::<f64, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        let jacobian = Jacobian2x2f::from(self.jacobian());
        let det_j = jacobian.determinant();

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

        let b_mat = h_mat * q_mat;
        // Gauss integration, area of standard tri is 0.5
        let core = b_mat.transpose() * elasticity_mat * b_mat * det_j;
        let stiffness_matrix: [[f64; 6]; 6] = (0.5 * t * core).into();
        stiffness_matrix
    }

    /// Get triangle element area value
    pub fn area(&self) -> f64 {
        let x = self.xs();
        let y = self.ys();
        let tri_area = 0.5 * ((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])).abs();
        tri_area
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
    fn k_printer(&mut self, material: (f64, f64)) {
        if self.k_matrix.is_none() {
            self.k_matrix = Some(self.calc_k(material));
        }
        print!("Tri2D3N k{} = \n[", self.id);
        for row in 0..6 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..6 {
                print!(" {:-9.6} ", self.k_matrix.unwrap()[row][col]);
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
