use crate::{node::Node2D, Dtype, K};
use na::*;
use std::fmt::{self, Write};

pub struct Rod2D2N<'rod> {
    pub id: usize,
    pub sec_area: Dtype,
    pub nodes: [&'rod Node2D; 2],
    pub k_matrix: Option<[[Dtype; 4]; 4]>,
}

impl<'rod> Rod2D2N<'rod> {
    /// Generate a 2D Rod2D2N element
    pub fn new(id: usize, sec_area: Dtype, nodes: [&Node2D; 2]) -> Rod2D2N {
        Rod2D2N {
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
            x_list[i] = self.nodes[i].coord[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in rod element
    pub fn ys(&self) -> [Dtype; 2] {
        let mut y_list = [0.0; 2];
        for i in 0..2 {
            y_list[i] = self.nodes[i].coord[1];
        }
        y_list
    }

    /// Get rod element length
    pub fn length(&self) -> (Dtype, Dtype, Dtype) {
        let dx = self.xs()[1] - self.xs()[0];
        let dy = self.ys()[1] - self.ys()[0];
        let l: Dtype = (dx * dx + dy * dy).sqrt();
        let cos: Dtype = ((dx as f64) / (l as f64)) as Dtype;
        let sin: Dtype = ((dy as f64) / (l as f64)) as Dtype;
        (l, cos, sin)
    }

    /// Get nodes' disps vector
    pub fn disps(&self) -> [Dtype; 4] {
        let mut disps = [0.0; 4];
        for idx in 0..2 {
            disps[2 * idx] = *self.nodes[idx].disps[0].borrow();
            disps[2 * idx + 1] = *self.nodes[idx].disps[1].borrow();
        }
        disps
    }

    /// Get nodes' force vector
    pub fn forces(&self) -> [Dtype; 4] {
        let mut forces = [0.0; 4];
        for idx in 0..2 {
            forces[2 * idx] = *self.nodes[idx].forces[0].borrow();
            forces[2 * idx + 1] = *self.nodes[idx].forces[1].borrow();
        }
        forces
    }

    /// Get transformation matrix
    pub fn trans_mat(&self) -> [[Dtype; 2]; 4] {
        let (_, cos, sin) = self.length();
        [[cos, 0.0], [sin, 0.0], [0.0, cos], [0.0, sin]]
    }

    /// Calculate element stiffness matrix K
    /// Return a 4x4 matrix, elements are Dtype
    fn calc_k(&self, material_args: (Dtype, Dtype)) -> [[Dtype; 4]; 4] {
        println!(
            "\n>>> Calculating Rod2D2N(#{})'s stiffness matrix k{}",
            self.id, self.id
        );
        let (ee, _nu) = material_args;
        let trans_mat = SMatrix::<Dtype, 2, 4>::from(self.trans_mat());
        let stiffness_mat: [[Dtype; 4]; 4] =
            ((ee * self.sec_area / self.length().0) * (trans_mat.transpose() * trans_mat)).into();
        stiffness_mat
    }

    /// Get element's strain vector, a scale in rod element
    fn calc_strain(&self) -> [Dtype; 3] {
        let unit = 1.0 / self.length().0;
        let b_mat = SMatrix::<Dtype, 1, 2>::from([-unit, unit]);
        let trans_mat = SMatrix::<Dtype, 2, 4>::from(self.trans_mat());
        let nodes_disps = SMatrix::<Dtype, 4, 1>::from(self.disps());
        let strain: [Dtype; 1] = (b_mat * trans_mat * nodes_disps).into();
        [strain[0], 0.0, 0.0]
    }

    /// Get element's strain vector, a scale in rod element
    fn calc_stress(&self, material_args: (Dtype, Dtype)) -> [Dtype; 3] {
        let (ee, _nu) = material_args;
        let stress: [Dtype; 1] = [ee * self.calc_strain()[0]];
        [stress[0], 0.0, 0.0]
    }
}

impl<'rod> K for Rod2D2N<'rod> {
    type Kmatrix = [[Dtype; 4]; 4];

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
                "!!! Rod2D2N#{}'s k mat is empty! call k() to calc it.",
                self.id
            );
        }

        print!("\nRod2D2N k{} =  (* 10^{})\n[", self.id, n_exp);
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
Element_1D Info:\n\tId:     {}\n\tArea:   {}\n\tLong:   {}\n\tType:   Rod2D2N\n\tNodes: {}\n\t       {}\n",
            self.id,
            self.sec_area,
            self.length().0,
            self.nodes[0],
            self.nodes[1]
        )
    }

    /// Get element's id number
    fn id(&self) -> usize {
        self.id
    }
}

impl fmt::Display for Rod2D2N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,"\n--------------------------------------------------------------------
Element_1D Info:\n\tId:     {}\n\tArea:   {}\n\tLong:   {}\n\tType:   Rod2D2N\n\tNodes: {}\n\t       {}\n",
self.id, self.sec_area, self.length().0, self.nodes[0], self.nodes[1])
    }
}
