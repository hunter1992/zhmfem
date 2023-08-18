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
}

impl fmt::Display for Rod2D2N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,"\n--------------------------------------------------------------------
Element_1D Info:\n\tId:     {}\n\tArea:   {}\n\tLong:   {}\n\tType:   Rod2D2N\n\tNodes: {}\n\t       {}\n",
self.id, self.sec_area, self.length().0, self.nodes[0], self.nodes[1])
    }
}
