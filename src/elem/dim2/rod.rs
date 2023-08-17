use crate::{node::Node2D, Dtype, K};
use na::*;
//use std::fmt::{self,Write};

pub struct Rod2D2N<'rod> {
    pub id: usize,
    pub alpha: Dtype,
    pub sec_area: Dtype,
    pub nodes: [&'rod Node2D; 2],
    pub k_matrix: Option<[[Dtype; 4]; 4]>,
}

impl<'rod> Rod2D2N<'rod> {
    /// Generate a 2D Rod2D2N element
    pub fn new(id: usize, alpha: Dtype, sec_area: Dtype, nodes: [&Node2D; 2]) -> Rod2D2N {
        Rod2D2N {
            id,
            alpha,
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
    pub fn length(&self) -> Dtype {
        let dx = self.xs()[0] - self.xs()[1];
        let dy = self.ys()[0] - self.ys()[1];
        (dx * dx + dy * dy).sqrt() as Dtype
    }

    /// Get nodes' disps vector
    pub fn disps(&self) -> [Dtype; 4] {
        let mut disps = [0.0; 4];
        for idx in 0..4 {
            disps[idx] = *self.nodes[idx].disps[0].borrow();
        }
        disps
    }

    /// Get nodes' force vector
    pub fn forces(&self) -> [Dtype; 4] {
        let mut forces = [0.0; 4];
        for idx in 0..4 {
            forces[idx] = *self.nodes[idx].forces[0].borrow();
        }
        forces
    }

    /// Get transformation matrix
    pub fn transform_to_global(&self) -> [[Dtype; 4]; 2] {
        let a: Dtype = self.alpha;
        let mut trans_mat: [[Dtype; 4]; 2] = [[0.0; 4]; 2];
        let cos = |x: Dtype| x.cos();
        let sin = |x: Dtype| x.sin();
        [[cos(a), sin(a), 0.0, 0.0], [0.0, 0.0, cos(a), sin(a)]]
    }
}
