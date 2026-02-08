use crate::dtty::{
    basic::{Dtype, Jacobian2D},
    matrix::CompressedMatrixSKS,
};
use crate::node::Node2D;
use crate::port::{SData, StaticStiffness};
use crate::tool::{compress_symmetry_matrix_sks, print_2darr};
use na::SMatrix;
use std::fmt::{self, Write};

/// Four-node triangular element in two-dimensional plane
/// The stress and strain in the quad element are linear
pub struct Quad2D4N<'quad2d4n> {
    pub id: usize,
    pub thick: Dtype,
    pub nodes: [&'quad2d4n Node2D; 4],
    pub k_matrix: Option<CompressedMatrixSKS>,
    pub material: [Dtype; 2],
}

impl<'quad2d4n> Quad2D4N<'quad2d4n> {
    /// Generate a new Quad2D4N element
    pub fn new(
        id: usize,
        thick: Dtype,
        nodes: [&'quad2d4n Node2D; 4],
        material: [Dtype; 2],
    ) -> Self {
        Quad2D4N {
            id,
            thick,
            nodes,
            k_matrix: None,
            material,
        }
    }

    /// Get element's id number
    pub fn id(&self) -> usize {
        self.id
    }

    /// Get the rectangle element area
    ///   D ----------- C
    ///     |         |
    ///     |         |
    ///     |         |
    ///   A ----------- B
    /// divide squdABCD into tow triangles: triABC and tri CDA
    /// S(ABCD) = S(ABC) + S(CDA)
    pub fn area(&self) -> Dtype {
        let x = self.get_nodes_xcoord();
        let y = self.get_nodes_ycoord();
        let dx1_21 = x[1] - x[0];
        let dx1_31 = x[2] - x[0];
        let dy1_21 = y[1] - y[0];
        let dy1_31 = y[2] - y[0];
        let dx2_21 = x[3] - x[2];
        let dx2_31 = x[0] - x[2];
        let dy2_21 = y[3] - y[2];
        let dy2_31 = y[0] - y[2];
        let s_tri1 = 0.5 * (dx1_21 * dy1_31 - dx1_31 * dy1_21).abs();
        let s_tri2 = 0.5 * (dx2_21 * dy2_31 - dx2_31 * dy2_21).abs();
        s_tri1 + s_tri2
    }

    /// Get the x-coords of nodes in quad element
    pub fn get_nodes_xcoord(&self) -> [Dtype; 4] {
        let mut x_list = [0.0; 4];
        for (idx, node) in self.nodes.iter().enumerate() {
            x_list[idx] = node.coords[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in quad element
    pub fn get_nodes_ycoord(&self) -> [Dtype; 4] {
        let mut y_list = [0.0; 4];
        for (idx, node) in self.nodes.iter().enumerate() {
            y_list[idx] = node.coords[1];
        }
        y_list
    }

    /// Get nodes' disps vector in quad element
    pub fn get_nodes_displacement(&self) -> [Dtype; 8] {
        let mut displacement = [0.0; 8];
        for (idx, node) in self.nodes.iter().enumerate() {
            displacement[2 * idx] = node.displs.borrow()[0];
            displacement[2 * idx + 1] = node.displs.borrow()[1];
        }
        displacement
    }

    /// Get nodes's force vector in quad element
    pub fn get_nodes_force(&self) -> [Dtype; 8] {
        let mut forces = [0.0; 8];
        for (idx, node) in self.nodes.iter().enumerate() {
            forces[2 * idx] = node.forces.borrow()[0];
            forces[2 * idx + 1] = node.forces.borrow()[1];
        }
        forces
    }

    /// Get any point's disps vector in quad element
    pub fn get_point_displacement(&self, s_t_coords: [Dtype; 2]) -> [Dtype; 2] {
        let s = s_t_coords[0];
        let t = s_t_coords[1];
        let n0 = self.shape_func_st(0usize)(s, t);
        let n1 = self.shape_func_st(1usize)(s, t);
        let n2 = self.shape_func_st(2usize)(s, t);
        let n3 = self.shape_func_st(3usize)(s, t);

        let disps = self.get_nodes_displacement();
        let u = n0 * disps[0] + n1 * disps[2] + n2 * disps[4] + n3 * disps[6];
        let v = n0 * disps[1] + n1 * disps[3] + n2 * disps[5] + n3 * disps[7];
        [u, v]
    }

    /// Get shape function matrix element N_i under Parametric coordinate
    /// The shape mat of quad elem:
    /// | N1   0    N2   0    N3   0    N4   0  |
    /// | 0    N1   0    N2   0    N3   0    N4 |
    ///
    /// N_i(xi, eta) = 0.25 * (1 + xi_i * xi) * (1 + eta_i * eta)
    ///
    /// -------------------
    ///   i  |  xi  |  eta
    /// -------------------
    ///   1  |  -1  |  -1   (从第三象限开始)
    ///   2  |   1  |  -1
    ///   3  |   1  |   1
    ///   4  |  -1  |   1
    /// -------------------
    ///
    /// N1 = 0.25 * (1 - s) * (1 - t)
    /// N2 = 0.25 * (1 + s) * (1 - t)
    /// N3 = 0.25 * (1 + s) * (1 + t)
    /// N4 = 0.25 * (1 - s) * (1 + t)
    fn shape_func_st(&self, ith: usize) -> impl Fn(Dtype, Dtype) -> Dtype {
        let s_sign = [-1.0, 1.0, 1.0, -1.0];
        let t_sign = [-1.0, -1.0, 1.0, 1.0];
        move |x: Dtype, y: Dtype| {
            0.25 * (1.0 + s_sign[ith] * x + t_sign[ith] * y + s_sign[ith] * t_sign[ith] * x * y)
        }
    }

    /// Calculate the Jacobian matrix of quadrilateral
    /// Jacobian = | dx/ds  dy/ds |
    ///            |              |
    ///            | dx/dt  dy/dt |
    fn jacobian(&self, s_t_coords: [Dtype; 2]) -> Jacobian2D {
        let x: [Dtype; 4] = self.get_nodes_xcoord();
        let y: [Dtype; 4] = self.get_nodes_ycoord();
        let xy = SMatrix::<Dtype, 4, 2>::from([x, y]);
        let dn_st = self.diff_shape_mat_st(s_t_coords);
        dn_st * xy
    }

    /// Derivatives of shape functions with respect to parametric coordinates
    /// 形函数对参数坐标的导数
    ///  | dN/ds |          | dN1/ds  dN2/ds  dN3/ds  dN4/ds |
    ///  |       |   =      |                                |
    ///  | dN/dt |          | dN1/dt  dN2/dt  dN3/dt  dN4/dt |
    ///
    ///                     | t-1     1-t     1+t     -1-t   |
    ///              = 0.25 |                                |
    ///                     | s-1    -1-s     1+s      1-s   |
    fn diff_shape_mat_st(&self, s_t_coords: [Dtype; 2]) -> SMatrix<Dtype, 2, 4> {
        let s: Dtype = s_t_coords[0];
        let t: Dtype = s_t_coords[1];
        0.25 * SMatrix::<Dtype, 4, 2>::from([
            [t - 1., -t + 1., t + 1., -t - 1.],
            [s - 1., -s - 1., s + 1., -s + 1.],
        ])
        .transpose()
    }

    /// Derivatives of shape functions with respect to physical coordinates
    /// 形函数对物理坐标的导数
    /// dN/dx = (dN/ds * ds/dx) + (dN/dt * dt/dx)
    /// dN/dy = (dN/ds * ds/dy) + (dN/dt * dt/dy)
    /// | dN/ds |        | dx/ds  dy/ds | | dN/dx |          | dN/dx |
    /// |       |    =   |              |*|       |   =  J * |       |
    /// | dN/dt |        | dx/dt  dy/dt | | dN/dy |2x4       | dN/dy |2x4
    ///
    /// | dN/dx |                 | dN/ds |
    /// |       |    =   J^(-1) * |       |
    /// | dN/dy |2x4              | dN/dt |2x4
    fn diff_shape_mat_xy(&self, s_t_coords: [Dtype; 2]) -> SMatrix<Dtype, 2, 4> {
        let diff_shape_mat_st = self.diff_shape_mat_st(s_t_coords);
        let jecobain = self.jacobian(s_t_coords);
        jecobain.try_inverse().unwrap() * diff_shape_mat_st
    }

    /// Get inverse shape matrix at integration points
    /// The shape mat at integration point is:
    ///     [ N1(s1, t1),   N2(s1, t1),   N3(s1, t1),   N4(s1, t1) ]
    ///     [ N1(s2, t2),   N2(s2, t2),   N3(s2, t2),   N4(s2, t2) ]
    ///     [ N1(s3, t3),   N2(s3, t3),   N3(s3, t3),   N4(s3, t3) ]
    ///     [ N1(s4, t4),   N2(s4, t4),   N3(s4, t4),   N4(s4, t4) ]
    pub fn inv_shape_mat_under_st_at_intpoint(&self) -> [[Dtype; 4]; 4] {
        let gp: Dtype = (1.0 as Dtype) / ((3.0 as Dtype).sqrt()); // gp for gauss point
        let n1 = self.shape_func_st(0);
        let n2 = self.shape_func_st(1);
        let n3 = self.shape_func_st(2);
        let n4 = self.shape_func_st(3);

        let shape_mat: SMatrix<Dtype, 4, 4> = SMatrix::<Dtype, 4, 4>::from([
            [n1(-gp, -gp), n1(gp, -gp), n1(gp, gp), n1(-gp, gp)],
            [n2(-gp, -gp), n2(gp, -gp), n2(gp, gp), n2(-gp, gp)],
            [n3(-gp, -gp), n3(gp, -gp), n3(gp, gp), n3(-gp, gp)],
            [n4(-gp, -gp), n4(gp, -gp), n4(gp, gp), n4(-gp, gp)],
        ]);
        shape_mat.try_inverse().unwrap().into()
    }

    /// Geometry matrix B(xi, eta) 用参数坐标表示物理坐标下的几何矩阵
    fn geometry_mat_xy(&self, s_t_coords: [Dtype; 2]) -> SMatrix<Dtype, 3, 8> {
        let dn_xy = self.diff_shape_mat_xy(s_t_coords);
        SMatrix::<Dtype, 3, 8>::from([
            [dn_xy[(0, 0)], 0.0, dn_xy[(1, 0)]],
            [0.0, dn_xy[(1, 0)], dn_xy[(0, 0)]],
            [dn_xy[(0, 1)], 0.0, dn_xy[(1, 1)]],
            [0.0, dn_xy[(1, 1)], dn_xy[(0, 1)]],
            [dn_xy[(0, 2)], 0.0, dn_xy[(1, 2)]],
            [0.0, dn_xy[(1, 2)], dn_xy[(0, 2)]],
            [dn_xy[(0, 3)], 0.0, dn_xy[(1, 3)]],
            [0.0, dn_xy[(1, 3)], dn_xy[(0, 3)]],
        ])
    }

    /// Calculate element stiffness matrix K
    /// return a 8x8 matrix, elements are Dtype
    pub fn calc_k(&self) -> [[Dtype; 8]; 8] {
        // println!(
        //     "\n>>> Calculating Quad2D4N(#{})'s stiffness matrix k{} ......",
        //     self.id, self.id
        // );
        let ee = self.material[0];
        let nu = self.material[1];
        let elasticity_mat = SMatrix::<Dtype, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        // 4 point Gauss integration, area of standard rectangle is 4.0
        let gauss_pt: Dtype = (1.0 as Dtype) / ((3.0 as Dtype).sqrt());
        let int_pts: [[[Dtype; 2]; 2]; 2] = [
            [[-gauss_pt, -gauss_pt], [gauss_pt, -gauss_pt]],
            [[gauss_pt, gauss_pt], [-gauss_pt, gauss_pt]],
        ];

        let mut k_matrix = SMatrix::<Dtype, 8, 8>::from([[0.0; 8]; 8]);
        // 2x2的4点Gauss积分，4个点权重都是1
        for row in 0..2 {
            for col in 0..2 {
                let j: Jacobian2D = self.jacobian(int_pts[row][col]);
                let det_j: Dtype = j.determinant();

                let b_mat = self.geometry_mat_xy(int_pts[row][col]);
                let core = b_mat.transpose() * elasticity_mat * b_mat * det_j;
                k_matrix += core;
            }
        }

        let stiffness_matrix: [[Dtype; 8]; 8] = (self.thick * k_matrix).into();
        stiffness_matrix
    }

    /// Get element's strain vector
    fn calc_strain(&self, s_t_coords: [Dtype; 2]) -> [Dtype; 3] {
        let b_mat = self.geometry_mat_xy(s_t_coords);
        let elem_nodes_disps = SMatrix::<Dtype, 8, 1>::from(self.get_nodes_displacement());
        let strain_vector: [Dtype; 3] = (b_mat * elem_nodes_disps).into();
        strain_vector
    }

    /// Get element's strss vector
    fn calc_stress(&self, s_t_coords: [Dtype; 2]) -> [Dtype; 3] {
        let ee = self.material[0];
        let nu = self.material[1];
        let elasticity_mat = SMatrix::<Dtype, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        let strain = SMatrix::<Dtype, 3, 1>::from(self.calc_strain(s_t_coords));
        let stress: [Dtype; 3] = (elasticity_mat * strain).into();
        stress
    }

    /// Get element integration points' strain
    pub fn calc_strain_at_intpts(&self) -> [[Dtype; 3]; 4] {
        let mut epsilon: [[Dtype; 3]; 4] = [[0.0; 3]; 4];
        let gauss_pt = (1.0 as Dtype) / ((3.0 as Dtype).sqrt());
        let int_pts: [[Dtype; 2]; 4] = [
            [-gauss_pt, -gauss_pt],
            [gauss_pt, -gauss_pt],
            [gauss_pt, gauss_pt],
            [-gauss_pt, gauss_pt],
        ];
        for idx in 0..4 {
            epsilon[idx] = self.calc_strain(int_pts[idx]);
        }
        epsilon
    }

    /// Get element integration points' strain
    pub fn calc_stress_at_intpts(&self) -> [[Dtype; 3]; 4] {
        let mut sigma: [[Dtype; 3]; 4] = [[0.0; 3]; 4];
        let gauss_pt = (1.0 as Dtype) / ((3.0 as Dtype).sqrt());
        let int_pts: [[Dtype; 2]; 4] = [
            [-gauss_pt, -gauss_pt],
            [gauss_pt, -gauss_pt],
            [gauss_pt, gauss_pt],
            [-gauss_pt, gauss_pt],
        ];
        for idx in 0..4 {
            sigma[idx] = self.calc_stress(int_pts[idx]);
        }
        sigma
    }

    /// Get element node's strain
    pub fn calc_strain_at_nodes(&self) -> [[Dtype; 3]; 4] {
        let s3: Dtype = (3.0 as Dtype).sqrt();
        let data: SMatrix<Dtype, 4, 3> =
            (SMatrix::<Dtype, 4, 4>::from([
                [1. + 0.5 * s3, -0.5, 1. - 0.5 * s3, -0.5],
                [-0.5, 1. + 0.5 * s3, -0.5, 1. - 0.5 * s3],
                [1. - 0.5 * s3, -0.5, 1. + 0.5 * s3, -0.5],
                [-0.5, 1. - 0.5 * s3, -0.5, 1. + 0.5 * s3],
            ])) * (SMatrix::<Dtype, 3, 4>::from(self.calc_strain_at_intpts()).transpose());

        let mut nodes_strain: [[Dtype; 3]; 4] = [[0.0; 3]; 4];
        for idx in 0..4 {
            for idy in 0..3 {
                nodes_strain[idx][idy] = data[(idx, idy)];
            }
        }
        nodes_strain
    }

    /// Get element node's stress
    pub fn calc_stress_at_nodes(&self) -> [[Dtype; 3]; 4] {
        let s3: Dtype = (3.0 as Dtype).sqrt();
        let data: SMatrix<Dtype, 4, 3> =
            (SMatrix::<Dtype, 4, 4>::from([
                [1. + 0.5 * s3, -0.5, 1. - 0.5 * s3, -0.5],
                [-0.5, 1. + 0.5 * s3, -0.5, 1. - 0.5 * s3],
                [1. - 0.5 * s3, -0.5, 1. + 0.5 * s3, -0.5],
                [-0.5, 1. - 0.5 * s3, -0.5, 1. + 0.5 * s3],
            ])) * (SMatrix::<Dtype, 3, 4>::from(self.calc_stress_at_intpts()).transpose());

        let mut nodes_stress: [[Dtype; 3]; 4] = [[0.0; 3]; 4];
        for idx in 0..4 {
            for idy in 0..3 {
                nodes_stress[idx][idy] = data[(idx, idy)];
            }
        }
        nodes_stress
    }

    /// Print element's strain value
    pub fn print_strain(&self, s_t_coords: [Dtype; 2]) {
        let strain = self.calc_strain(s_t_coords);
        println!(
            "\nelem[{}] strain:\n\tE_xx = {:-16.6}\n\tE_yy = {:-16.6}\n\tE_xy = {:-16.6}",
            self.id, strain[0], strain[1], strain[2]
        );
    }

    /// Print element's stress value
    pub fn print_stress(&self, s_t_coords: [Dtype; 2]) {
        let stress = self.calc_stress(s_t_coords);
        println!(
            "\nelem[{}] stress:\n\tS_xx = {:-16.6}\n\tS_yy = {:-16.6}\n\tS_xy = {:-16.6}",
            self.id, stress[0], stress[1], stress[2]
        );
    }
}

impl<'quad2d4n> StaticStiffness for Quad2D4N<'quad2d4n> {
    /// Cache stiffness matrix for quad element
    fn k(&mut self) -> &CompressedMatrixSKS {
        if self.k_matrix.is_none() {
            self.k_matrix
                .get_or_insert(compress_symmetry_matrix_sks(&self.calc_k()))
        } else {
            self.k_matrix.as_ref().unwrap()
        }
    }

    /// Print quad element's stiffness matrix
    fn k_printr(&self, n_exp: Dtype) {
        print_2darr(
            "\n Tri2D4N k",
            self.id,
            &self.k_matrix.as_ref().unwrap().recover_square_arr::<8>(),
            n_exp,
        );
    }

    /// Return quad elem's stiffness matrix's format string
    fn k_string(&self, n_exp: Dtype) -> String {
        let mut k_matrix = String::new();
        let elem_stiffness_mat: [[Dtype; 8]; 8] =
            self.k_matrix.as_ref().unwrap().recover_square_arr::<8>();
        for row in 0..8 {
            if row == 0 {
                write!(k_matrix, "[[").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, " [").expect("!!! Write tri k_mat failed!");
            }
            for col in 0..8 {
                write!(
                    k_matrix,
                    " {:>-12.6} ",
                    elem_stiffness_mat[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                )
                .expect("!!! Write tri k_mat failed!");
            }
            if row == 7 {
                write!(k_matrix, "]]").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, "]\n").expect("!!! Write tri k_mat failed!");
            }
        }
        k_matrix
    }
}

impl<'quad2d4n> SData for Quad2D4N<'quad2d4n> {
    fn elem_size(&self) -> Dtype {
        self.area()
    }

    fn nodes_ids(&self) -> Vec<usize> {
        let mut id: Vec<usize> = Vec::with_capacity(4);
        for idx in 0..4 {
            id.push(self.nodes[idx].id);
        }
        id
    }

    fn strain_at_nodes(&mut self) -> Vec<Dtype> {
        let node_strain = self.calc_strain_at_nodes();
        node_strain
            .into_iter()
            .flat_map(|value| value.into_iter())
            .collect()
    }

    fn stress_at_nodes(&mut self) -> Vec<Dtype> {
        let node_stress = self.calc_stress_at_nodes();
        node_stress
            .into_iter()
            .flat_map(|value| value.into_iter())
            .collect()
    }
}

impl fmt::Display for Quad2D4N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let gauss_pt = (1.0 as Dtype) / ((3.0 as Dtype).sqrt());
        let int_pts: [[[Dtype; 2]; 2]; 2] = [
            [[-gauss_pt, -gauss_pt], [gauss_pt, -gauss_pt]],
            [[gauss_pt, gauss_pt], [-gauss_pt, gauss_pt]],
        ];
        write!(
            f,
            "\n-----------------------------------------------------------------------------\nElem_Quad2D4N:\n\tId:\t{}\n\tArea: {:-12.6}\n\tMats: {:-12.6} (Young's modulus)\n\t      {:-12.6} (Poisson's ratio)\n\tNodes:{}{}{}{}\n\tStrain: (at integration points)\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\n\tStress: (at integration points)\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\n\tStiffness Matrix K{} =  (*10^0)\n{}",
            self.id,
            self.area(),
            self.material[0],
            self.material[1],
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
            self.nodes[3],
            self.calc_strain(int_pts[0][0]),
            self.calc_strain(int_pts[0][1]),
            self.calc_strain(int_pts[1][0]),
            self.calc_strain(int_pts[1][1]),
            self.calc_stress(int_pts[0][0]),
            self.calc_stress(int_pts[0][1]),
            self.calc_stress(int_pts[1][0]),
            self.calc_stress(int_pts[1][1]),
            self.id,
            self.k_string(0.0)
        )
    }
}
