use crate::dtty::{
    basic::{Dtype, Jacobian3D},
    matrix::CompressedMatrixSKS,
};
use crate::node::Node3D;
use crate::port::{SData, StaticStiffness};
use crate::tool::{compress_symmetry_matrix_sks, print_2darr};
use na::SMatrix;
use std::arch::x86_64::*;
use std::fmt::{self, Write};

/// Eight-node hexahedron element in 3D space
pub struct Hex3D8N<'hex3d8n> {
    pub id: usize,
    pub nodes: [&'hex3d8n Node3D; 8],
    pub k_matrix: Option<CompressedMatrixSKS>,
    pub material: [Dtype; 2],
}

impl<'hex3d8n> Hex3D8N<'hex3d8n> {
    /// Generate a new Quad2D4N element
    pub fn new(id: usize, nodes: [&'hex3d8n Node3D; 8], material: [Dtype; 2]) -> Self {
        Hex3D8N {
            id,
            nodes,
            k_matrix: None,
            material,
        }
    }

    /// Get element's id number
    pub fn id(&self) -> usize {
        self.id
    }

    /// Get the x-coords of nodes in quad element
    pub fn get_nodes_xcoord(&self) -> [Dtype; 8] {
        let mut x_list = [0.0; 8];
        for (idx, node) in self.nodes.iter().enumerate() {
            x_list[idx] = node.coords[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in quad element
    pub fn get_nodes_ycoord(&self) -> [Dtype; 8] {
        let mut y_list = [0.0; 8];
        for (idx, node) in self.nodes.iter().enumerate() {
            y_list[idx] = node.coords[1];
        }
        y_list
    }

    /// Get the z-coords of nodes in quad element
    pub fn get_nodes_zcoord(&self) -> [Dtype; 8] {
        let mut z_list = [0.0; 8];
        for (idx, node) in self.nodes.iter().enumerate() {
            z_list[idx] = node.coords[2];
        }
        z_list
    }

    /// Get nodes' disps vector in quad element
    pub fn get_nodes_displacement(&self) -> [Dtype; 24] {
        let mut displacement = [0.0; 24];
        for (idx, node) in self.nodes.iter().enumerate() {
            displacement[3 * idx] = node.displs.borrow()[0];
            displacement[3 * idx + 1] = node.displs.borrow()[1];
            displacement[3 * idx + 2] = node.displs.borrow()[2];
        }
        displacement
    }

    /// Get nodes's force vector in quad element
    pub fn get_nodes_force(&self) -> [Dtype; 24] {
        let mut forces = [0.0; 24];
        for (idx, node) in self.nodes.iter().enumerate() {
            forces[3 * idx] = node.forces.borrow()[0];
            forces[3 * idx + 1] = node.forces.borrow()[1];
            forces[3 * idx + 2] = node.forces.borrow()[2];
        }
        forces
    }

    /// Hex3D8N element shape function matrix under parametric coords
    /// The shape mat of quad elem:
    /// | N1  0  0  N2  0  0  N3  0  0  N4  0  0  N5  0  0  N6  0  0  N7  0  0  N8  0  0 |
    /// | 0  N1  0  0  N2  0  0  N3  0  0  N4  0  0  N5  0  0  N6  0  0  N7  0  0  N8  0 |
    /// | 0   0 N1  0  0  N2  0  0  N3  0  0  N4  0  0  N5  0  0  N6  0  0  N7  0  0  N8 |
    ///
    /// N_i(xi, eta, zeta) = 0.125 * (1 + xi_i * xi)
    ///                            * (1 + eta_i * eta)          (i = 1,...,8)
    ///                            * (1 + zeta_i * zeta)
    ///
    /// ----------------------------
    ///   i  |  xi  |  eta  |  zeta
    /// ----------------------------
    ///   1  |  -1  |  -1   |   -1      (从下层第四象限开始)
    ///   2  |   1  |  -1   |   -1
    ///   3  |   1  |   1   |   -1
    ///   4  |  -1  |   1   |   -1
    ///   5  |  -1  |  -1   |    1
    ///   6  |   1  |  -1   |    1
    ///   7  |   1  |   1   |    1
    ///   8  |  -1  |   1   |    1
    /// ----------------------------
    ///
    /// N1 = 0.125 * (1 - r) * (1 - s) * (1 - t)
    /// N2 = 0.125 * (1 + r) * (1 - s) * (1 - t)
    /// N3 = 0.125 * (1 + r) * (1 + s) * (1 - t)
    /// N4 = 0.125 * (1 - r) * (1 + s) * (1 - t)
    /// N5 = 0.125 * (1 - r) * (1 - s) * (1 + t)
    /// N6 = 0.125 * (1 + r) * (1 - s) * (1 + t)
    /// N7 = 0.125 * (1 + r) * (1 + s) * (1 + t)
    /// N8 = 0.125 * (1 - r) * (1 + s) * (1 + t)
    fn shape_func_rst(&self, ith: usize) -> impl Fn(Dtype, Dtype, Dtype) -> Dtype {
        let r_sign: [Dtype; 8] = [-1., 1., 1., -1., -1., 1., 1., -1.];
        let s_sign: [Dtype; 8] = [-1., -1., 1., 1., -1., -1., 1., 1.];
        let t_sign: [Dtype; 8] = [-1., -1., -1., -1., 1., 1., 1., 1.];
        move |x: Dtype, y: Dtype, z: Dtype| {
            0.125 * (1.0 + r_sign[ith] * x) * (1.0 + s_sign[ith] * y) * (1.0 + t_sign[ith] * z)
        }
    }
}

impl fmt::Display for Hex3D8N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\n-----------------------------------------------------------------------------"
        )
    }
}
