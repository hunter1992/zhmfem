use crate::{node::Node1D, K};
use na::*;
use std::fmt::{self, Write};

pub struct Rod1D2N<'rod> {
    pub id: usize,
    pub sec_area: f64,
    pub nodes: [&'rod Node1D; 2],
    pub k_matrix: Option<[[f64; 2]; 2]>,
}

impl<'rod> Rod1D2N<'rod> {
    /// Generate a 1D Rod1D2N element
    pub fn new(id: usize, sec_area: f64, nodes: [&Node1D; 2]) -> Rod1D2N {
        Rod1D2N {
            id,
            sec_area,
            nodes,
            k_matrix: None,
        }
    }

    /// Get the x-coords of nodes in rod element
    pub fn xs(&self) -> [f64; 2] {
        let mut x_list = [0.0; 2];
        for i in 0..2 {
            x_list[i] = self.nodes[i].coord[0];
        }
        x_list
    }

    /// Get rod element length
    pub fn length(&self) -> f64 {
        let x = self.xs();
        (x[0] - x[1]).abs()
    }

    /// Get nodes' disps vector
    pub fn disps(&self) -> [f64; 2] {
        let mut disps = [0.0; 2];
        for idx in 0..2 {
            disps[idx] = *self.nodes[idx].disps[0].borrow();
        }
        disps
    }

    /// Get nodes' force vector
    pub fn forces(&self) -> [f64; 2] {
        let mut forces = [0.0; 2];
        for idx in 0..2 {
            forces[idx] = *self.nodes[idx].forces[0].borrow();
        }
        forces
    }

    /// Get point's disp in element coord
    pub fn point_disp(&self, point_coord: [f64; 1]) -> [f64; 1] {
        let x = point_coord[0];
        let n0 = self.shape_mat_i(0usize)(x);
        let n1 = self.shape_mat_i(1usize)(x);

        let node_disps = self.disps();
        let u = n0 * node_disps[0] + n1 * node_disps[1];
        [u]
    }

    /// Get shape matrix element N_i
    fn shape_mat_i(&self, i: usize) -> impl Fn(f64) -> f64 {
        /* The shape mat of rod elem:
         * [N1 N2]  which is a 1x2 mat */
        let a: [f64; 2] = [1.0, 0.0];
        let b: [f64; 2] = [-1.0, 1.0];
        let length = self.length();

        move |x: f64| a[i] + b[i] * x / length
    }

    /// Calculate element stiffness matrix K
    /// Return a 2x2 matrix, elements are f64
    fn calc_k(&self, material_args: (f64, f64)) -> [[f64; 2]; 2] {
        println!(
            "\n>>> Calculating Rod1D2N(#{})'s stiffness matrix k{} ......",
            self.id, self.id
        );
        let (ee, _nu) = material_args;
        let stiffness_matrix: [[f64; 2]; 2] =
            (SMatrix::<f64, 2, 2>::from([[1.0, -1.0], [-1.0, 1.0]])
                * (ee * self.sec_area / self.length()))
            .into();
        stiffness_matrix
    }

    /// Get element's strain vector, in 1d it's a scale
    pub fn strain(&self) -> [f64; 1] {
        let iden: f64 = 1.0 / self.length();
        let node_disps = self.disps();
        let strain: [f64; 1] = [-iden * node_disps[0] + iden * node_disps[1]];
        strain
    }

    /// Get element's stress vector, in 1d it's a scale
    pub fn stress(&self, material_args: (f64, f64)) -> [f64; 1] {
        let (ee, _nu) = material_args;
        let stress: [f64; 1] = [ee * self.strain()[0]];
        stress
    }

    /// Get element's info string
    pub fn info(&self) -> String {
        format!(
            "\n--------------------------------------------------------------------
Element_1D Info:\n\tId:     {}\n\tArea:   {}\n\tLong:   {}\n\tType:   Rod1D2N\n\tNodes: {}\n\t       {}\n",
            self.id,
            self.sec_area,
            self.length(),
            self.nodes[0],
            self.nodes[1]
        )
    }
}

/// Implement zhm::K trait for triangle element
impl<'rod> K for Rod1D2N<'rod> {
    type Kmatrix = [[f64; 2]; 2];

    /// Cache stiffness matrix for rod element
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

    /// Print rod element's stiffness matrix
    fn k_printer(&self, n_exp: f64) {
        if self.k_matrix.is_none() {
            panic!(
                "!!! Rod1D2N#{}'s k mat is empty! call k() to calc it.",
                self.id
            );
        }

        print!("\nRod1D2N k{} =  (* 10^{})\n[", self.id, n_exp as u8);
        for row in 0..2 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..2 {
                print!(
                    " {:>-10.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp))
                );
            }
            if row == 1 {
                println!("]]");
            } else {
                println!("]");
            }
        }
        print!("\n");
    }

    /// Return rod elem's stiffness matrix's format string
    fn k_string(&self, n_exp: f64) -> String {
        let mut k_matrix = String::new();
        for row in 0..2 {
            if row == 0 {
                write!(k_matrix, "[[").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, " [").expect("!!! Write tri k_mat failed!");
            }
            for col in 0..2 {
                write!(
                    k_matrix,
                    " {:>-10.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp))
                )
                .expect("!!! Write tri k_mat failed!");
            }
            if row == 1 {
                write!(k_matrix, "]]").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, "]\n").expect("!!! Write tri k_mat failed!");
            }
        }
        k_matrix
    }
}

impl fmt::Display for Rod1D2N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "\n--------------------------------------------------------------------
Element_1D Info:\n\tId:     {}\n\tArea:   {}\n\tLong:   {}\n\tType:   Rod1D2N\n\tNodes: {}\n\t       {}\n",
              self.id,
              self.sec_area,
              self.length(),
              self.nodes[0],
              self.nodes[1])
    }
}
