extern crate nalgebra as na;

use super::triangle::Tri2D3N;
use crate::{node::Node2D, K};
use na::*;
use std::fmt;

pub struct Quad2D4N<'quad> {
    pub id: usize,
    pub thick: f64,
    pub nodes: [&'quad Node2D; 4],
    pub k_matrix: Option<[[f64; 8]; 8]>,
}

impl<'quad> Quad2D4N<'quad> {
    /// generate a new Rec2D4N element
    pub fn new(id: usize, thick: f64, nodes: [&Node2D; 4]) -> Quad2D4N {
        Quad2D4N {
            id,
            thick,
            nodes,
            k_matrix: None,
        }
    }

    /// get the x coords of nodes in Rec2D4N element
    pub fn xs(&self) -> [f64; 4] {
        let mut x_list = [0.0; 4];
        for i in 0..4 {
            x_list[i] = self.nodes[i].coord[0];
        }
        x_list
    }

    /// get the y coords of nodes in tri element
    pub fn ys(&self) -> [f64; 4] {
        let mut y_list = [0.0; 4];
        for i in 0..4 {
            y_list[i] = self.nodes[i].coord[1];
        }
        y_list
    }

    /// Calculate the Jacobian matrix of quadrilateral
    pub fn jacobian(&self, xi_eta: [f64; 2]) -> [[f64; 2]; 2] {
        let x: [f64; 4] = self.xs();
        let y: [f64; 4] = self.ys();
        let xi: f64 = xi_eta[0];
        let eta: f64 = xi_eta[1];

        let p_xi: f64 = 1.0 + xi;
        let n_xi: f64 = 1.0 - xi;
        let p_eta: f64 = 1.0 + eta;
        let n_eta: f64 = 1.0 - eta;

        let x_d_xi: f64 = 0.25 * (-n_eta * x[0] + n_eta * x[1] + p_eta * x[2] - p_eta * x[3]);
        let y_d_xi: f64 = 0.25 * (-n_eta * y[0] + n_eta * y[1] + p_eta * y[2] - p_eta * y[3]);
        let x_d_eta: f64 = 0.25 * (-n_xi * x[0] - p_xi * x[1] + p_xi * x[2] + n_xi * x[3]);
        let y_d_eta: f64 = 0.25 * (-n_xi * y[0] - p_xi * y[1] + p_xi * y[2] + n_xi * y[3]);
        [[x_d_xi, y_d_xi], [x_d_eta, y_d_eta]]
    }

    /// Calculate element stiffness matrix K
    /// return a 8x8 matrix, elements are f64
    fn calc_k(&self, material_args: (f64, f64)) -> [[f64; 8]; 8] {
        println!(
            "\n>>> Calculating Quad2D4N#{}'s stiffness matrix k{} ......",
            self.id, self.id
        );
        let (ee, nu) = material_args;
        let elasticity_mat = SMatrix::<f64, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        let h_mat = |j: [[f64; 2]; 2], det_j: f64| {
            SMatrix::<f64, 3, 4>::from([
                [j[1][1], 0.0, -j[1][0]],
                [-j[0][1], 0.0, j[0][0]],
                [0.0, -j[1][0], j[1][1]],
                [0.0, j[0][0], -j[0][1]],
            ]) / det_j
        };

        let q_mat = |xi_eta: [f64; 2]| {
            let p_xi: f64 = 1.0 + xi_eta[0];
            let n_xi: f64 = 1.0 - xi_eta[0];
            let p_eta: f64 = 1.0 + xi_eta[1];
            let n_eta: f64 = 1.0 - xi_eta[1];

            SMatrix::<f64, 4, 8>::from([
                [-n_eta, -n_xi, 0.0, 0.0],
                [0.0, 0.0, -n_eta, -n_xi],
                [n_eta, -p_xi, 0.0, 0.0],
                [0.0, 0.0, n_eta, -p_xi],
                [p_eta, p_xi, 0.0, 0.0],
                [0.0, 0.0, p_eta, p_xi],
                [-p_eta, n_xi, 0.0, 0.0],
                [0.0, 0.0, -p_eta, n_xi],
            ]) * 0.25
        };

        // 4 point Gauss integration, area of standard rectangle is 4.0
        let gauss_pt = 1.0f64 / (3.0f64.sqrt());
        let int_pts = [
            [[-gauss_pt, -gauss_pt], [gauss_pt, -gauss_pt]],
            [[gauss_pt, gauss_pt], [-gauss_pt, gauss_pt]],
        ];

        let mut k_matrix = SMatrix::<f64, 8, 8>::from([[0.0; 8]; 8]);

        for row in 0..2 {
            for col in 0..2 {
                let j_raw = self.jacobian(int_pts[row][col]);
                let j = SMatrix::<f64, 2, 2>::from(j_raw);
                let det_j = j.determinant().abs();

                let b_mat = h_mat(j_raw, det_j) * q_mat(int_pts[row][col]);
                let core = b_mat.transpose() * elasticity_mat * b_mat * det_j;
                k_matrix += core;
            }
        }

        let stiffness_matrix: [[f64; 8]; 8] = (self.thick * k_matrix).into();
        stiffness_matrix
    }

    /// Get the rectangle element area
    pub fn area(&self) -> f64 {
        let tri1: Tri2D3N = Tri2D3N {
            id: 0,
            thick: self.thick,
            nodes: [self.nodes[0], self.nodes[1], self.nodes[2]],
            k_matrix: None,
        };
        let tri2: Tri2D3N = Tri2D3N {
            id: 0,
            thick: self.thick,
            nodes: [self.nodes[3], self.nodes[1], self.nodes[2]],
            k_matrix: None,
        };
        let rect_area = tri1.area() + tri2.area();
        rect_area
    }
}

impl<'quad> K for Quad2D4N<'quad> {
    type Kmatrix = [[f64; 8]; 8];

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

    fn k_printer(&mut self, material: (f64, f64)) {
        if self.k_matrix.is_none() {
            self.k_matrix = Some(self.calc_k(material));
        }
        print!("Quad2D4N k{} = \n[", self.id);
        for row in 0..8 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..8 {
                print!(" {:-9.6} ", self.k_matrix.unwrap()[row][col]);
            }
            if row == 7 {
                println!("]]");
            } else {
                println!("]");
            }
        }
        println!("");
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
