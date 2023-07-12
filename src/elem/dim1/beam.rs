use crate::{node::Node2D, Dtype, K};
use na::*;
use std::fmt::{self, Write};

/// 细长梁以及纯弯梁，抽象成具有如下特征的理论模型：
///   1) 仅用x坐标刻画所有位置(即细长之含义)；
///   2) 变形主要为垂直于x轴的挠度
/// 但是，节点的自由度不为1，而包含：挠度(y)、转角theta，
/// 此亦是梁单元为非协调单元的原因。
pub struct Beam1D2N<'beam1d2n> {
    pub id: usize,
    pub sec_area: Dtype,
    pub nodes: [&'beam1d2n Node2D; 2],
    pub k_matrix: Option<[[Dtype; 4]; 4]>,
}

impl<'beam1d2n> Beam1D2N<'beam1d2n> {
    /// Generate a Beam1D2N element
    pub fn new(id: usize, sec_area: Dtype, nodes: [&Node2D; 2]) -> Beam1D2N {
        Beam1D2N {
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

    /// Get the x-coords of nodes in Beam1D2N element
    pub fn xs(&self) -> [Dtype; 2] {
        let mut x_list = [0.0; 2];
        for i in 0..2 {
            x_list[i] = self.nodes[i].coord[0];
        }
        x_list
    }

    /// Get rod element length
    pub fn length(&self) -> Dtype {
        let x = self.xs();
        (x[0] - x[1]).abs()
    }

    /// Get nodes' displacement/rotate angle vector
    pub fn disps(&self) -> [Dtype; 4] {
        let mut disps = [0.0; 4];
        for idx in 0..2 {
            disps[2 * idx] = *self.nodes[idx].disps[0].borrow();
            disps[2 * idx + 1] = *self.nodes[idx].disps[1].borrow();
        }
        disps
    }

    /// Get nodes' force/moment vector
    pub fn forces(&self) -> [Dtype; 4] {
        let mut forces = [0.0; 4];
        for idx in 0..2 {
            forces[2 * idx] = *self.nodes[idx].forces[0].borrow();
            forces[2 * idx + 1] = *self.nodes[idx].forces[1].borrow();
        }
        forces
    }
}
