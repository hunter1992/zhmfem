use crate::{node::Node1D, Dtype, K};
use na::*;
use std::fmt::{self, Write};

pub struct Rod1D2N<'rod> {
    pub id: usize,
    pub sec_area: Dtype,
    pub nodes: [&'rod Node1D; 2],
    pub k_matrix: Option<[[Dtype; 2]; 2]>,
}

impl<'rod> Rod1D2N<'rod> {
    /// Generate a 1D Rod1D2N element
    pub fn new(id: usize, sec_area: Dtype, nodes: [&Node1D; 2]) -> Rod1D2N {
        Rod1D2N {
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

    /// Get the x-coords of nodes in rod element
    pub fn xs(&self) -> [Dtype; 2] {
        let mut x_list = [0.0; 2];
        for i in 0..2 {
            x_list[i] = self.nodes[i].coords[0];
        }
        x_list
    }

    /// Get rod element length
    pub fn length(&self) -> Dtype {
        let x = self.xs();
        (x[0] - x[1]).abs()
    }

    /// Get nodes' disps vector
    pub fn disps(&self) -> [Dtype; 2] {
        let mut disps = [0.0; 2];
        for idx in 0..2 {
            disps[idx] = self.nodes[idx].displs[0];
        }
        disps
    }

    /// Get nodes' force vector
    pub fn forces(&self) -> [Dtype; 2] {
        let mut forces = [0.0; 2];
        for idx in 0..2 {
            forces[idx] = self.nodes[idx].forces[0];
        }
        forces
    }

    /// Get point's disp in element coord
    pub fn point_disp(&self, point_coord: [Dtype; 1]) -> [Dtype; 1] {
        let x = point_coord[0];
        let n0 = self.shape_mat_i(0usize)(x);
        let n1 = self.shape_mat_i(1usize)(x);

        let node_disps = self.disps();
        let u = n0 * node_disps[0] + n1 * node_disps[1];
        [u]
    }

    /// Get shape matrix element N_i
    fn shape_mat_i(&self, i: usize) -> impl Fn(Dtype) -> Dtype {
        /* The shape mat of rod elem:
         * [N1 N2]  which is a 1x2 mat
         * N1 = 1 - x/L = 1 - epsilon
         * N2 = x/L     = epsilon      */
        let a: [Dtype; 2] = [1.0, 0.0];
        let b: [Dtype; 2] = [-1.0, 1.0];
        let length = self.length();
        // a[0] + b[0]*epsilon构造出N1,a[1] + b[1]*epsilon构造出N2
        move |x: Dtype| a[i] + b[i] * x / length
    }

    /// Calculate element stiffness matrix K
    /// Return a 2x2 matrix, elements are Dtype
    fn calc_k(&self, material_args: (Dtype, Dtype)) -> [[Dtype; 2]; 2] {
        println!(
            "\n>>> Calculating Rod1D2N(#{})'s stiffness matrix k{} ......",
            self.id, self.id
        );
        let (ee, _nu) = material_args;
        let stiffness_matrix: [[Dtype; 2]; 2] =
            (SMatrix::<Dtype, 2, 2>::from([[1.0, -1.0], [-1.0, 1.0]])
                * (ee * self.sec_area / self.length()))
            .into();
        stiffness_matrix
    }

    /// Get element's strain vector, in 1d it's a scale
    fn calc_strain(&self) -> [Dtype; 1] {
        let unit: Dtype = 1.0 / self.length();
        let node_disps = self.disps();
        let strain: [Dtype; 1] = [-unit * node_disps[0] + unit * node_disps[1]];
        strain
    }

    /// Get element's stress vector, in 1d it's a scale
    fn calc_stress(&self, material_args: (Dtype, Dtype)) -> [Dtype; 1] {
        let (ee, _nu) = material_args;
        let stress: [Dtype; 1] = [ee * self.calc_strain()[0]];
        stress
    }
}

/// Implement zhm::K trait for triangle element
impl<'rod> K for Rod1D2N<'rod> {
    type Kmatrix = [[Dtype; 2]; 2];

    /// Cache stiffness matrix for rod element
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

    /// Print rod element's stiffness matrix
    fn k_printer(&self, n_exp: Dtype) {
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
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
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
    fn k_string(&self, n_exp: Dtype) -> String {
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
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
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

    /// Get the stress at (x) inside the element
    fn strain(&self, _xyz: [Dtype; 3]) -> Vec<Dtype> {
        self.calc_strain().to_vec()
    }

    /// Get the stress at (x) inside the element
    fn stress(&self, _xyz: [Dtype; 3], material: (Dtype, Dtype)) -> Vec<Dtype> {
        self.calc_stress(material).to_vec()
    }

    /// Get element's info string
    fn info(&self) -> String {
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

    /// Get element's id number
    fn id(&self) -> usize {
        self.id
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
