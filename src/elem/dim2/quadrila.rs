use super::triangle::Tri2D3N;
use crate::dtty::{
    basic::{Dtype, Jacobian2D},
    matrix::CompressedMatrix,
};
use crate::node::Node2D;
use crate::port::K;
use crate::tool::{compress_matrix, print_1darr};
use na::{SMatrix, SVector};
use std::fmt::{self, Write};

pub struct Quad2D4N<'quad2d4n> {
    pub id: usize,
    pub thick: Dtype,
    pub nodes: [&'quad2d4n Node2D; 4],
    pub strain: Option<[[Dtype; 3]; 4]>,
    pub stress: Option<[[Dtype; 3]; 4]>,
    pub k_matrix: Option<CompressedMatrix>,
    pub material: &'quad2d4n (Dtype, Dtype),
}

impl<'quad2d4n> Quad2D4N<'quad2d4n> {
    /// Generate a new Quad2D4N element
    pub fn new(
        id: usize,
        thick: Dtype,
        nodes: [&'quad2d4n Node2D; 4],
        material: &'quad2d4n (Dtype, Dtype),
    ) -> Self {
        Quad2D4N {
            id,
            thick,
            nodes,
            strain: None,
            stress: None,
            k_matrix: None,
            material,
        }
    }

    /// Set element material_args
    pub fn set_material(&mut self, material_args: &'quad2d4n (Dtype, Dtype)) {
        self.material = material_args;
    }

    /// Get the rectangle element area
    pub fn area(&self) -> Dtype {
        let tri1: Tri2D3N = Tri2D3N {
            id: 0,
            thick: self.thick,
            nodes: [self.nodes[0], self.nodes[1], self.nodes[2]],
            strain: None,
            stress: None,
            k_matrix: None,
            material: self.material,
        };
        let tri2: Tri2D3N = Tri2D3N {
            id: 1,
            thick: self.thick,
            nodes: [self.nodes[2], self.nodes[3], self.nodes[0]],
            strain: None,
            stress: None,
            k_matrix: None,
            material: self.material,
        };
        tri1.area() + tri2.area()
    }

    /// Get the x-coords of nodes in tri element
    pub fn get_nodes_xcoords(&self) -> [Dtype; 4] {
        let mut x_list = [0.0; 4];
        for i in 0..4 {
            x_list[i] = self.nodes[i].coords[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in tri element
    pub fn get_nodes_ycoords(&self) -> [Dtype; 4] {
        let mut y_list = [0.0; 4];
        for i in 0..4 {
            y_list[i] = self.nodes[i].coords[1];
        }
        y_list
    }

    /// Get nodes' disps vector in tri element
    pub fn get_nodes_displacement(&self) -> [Dtype; 8] {
        let mut disps = [0.0; 8];
        for idx in 0..4 {
            disps[2 * idx] = self.nodes[idx].displs.borrow()[0];
            disps[2 * idx + 1] = self.nodes[idx].displs.borrow()[1];
        }
        disps
    }

    /// Get nodes's force vector in tri element
    pub fn get_nodes_force(&self) -> [Dtype; 8] {
        let mut forces = [0.0; 8];
        for idx in 0..4 {
            forces[2 * idx] = self.nodes[idx].forces.borrow()[0];
            forces[2 * idx + 1] = self.nodes[idx].forces.borrow()[1];
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

    /// Get shape function matrix element N_i
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
    /// N_i(xi, eta) = 0.25 * (1 + xi_i * xi) * (1 + eta_i * eta)
    fn shape_func_st(&self, ith: usize) -> impl Fn(Dtype, Dtype) -> Dtype {
        /* The shape mat of quad elem:
         * |N1   0    N2   0    N3   0    N4   0 |
         * |0    N1   0    N2   0    N3   0    N4|
         *
         * N1 = 0.25 * (1 - Xi * s) * (1 - Yi * t)
         * N2 = 0.25 * (1 + Xi * s) * (1 - Yi * t)
         * N3 = 0.25 * (1 + Xi * s) * (1 + Yi * t)
         * N4 = 0.25 * (1 - Xi * s) * (1 + Yi * t)
         *
         */
        let s_sign = [-1.0, 1.0, 1.0, -1.0];
        let t_sign = [-1.0, -1.0, 1.0, 1.0];
        move |x: Dtype, y: Dtype| {
            0.25 * (1.0 + s_sign[ith] * x + t_sign[ith] * y + s_sign[ith] * t_sign[ith] * x * y)
        }
    }

    /// Derivatives of shape functions with respect to parametric coordinates
    /// 形函数对参数坐标的导数
    ///  | dN/ds |          | dN1/ds  dN1/ds  dN1/ds  dN1/ds |
    ///  |       |   =      |                                |
    ///  | dN/dt |          | dN1/dt  dN1/dt  dN1/dt  dN1/dt |
    ///
    ///                     | t-1     1-t     1+t     -1-t   |
    ///              = 0.25 |                                |
    ///                     | s-1     -1-s    1+s     1-s    |
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

    /// Calculate the Jacobian matrix of quadrilateral
    fn jacobian(&self, s_t_coords: [Dtype; 2]) -> Jacobian2D {
        let x: [Dtype; 4] = self.get_nodes_xcoords();
        let y: [Dtype; 4] = self.get_nodes_ycoords();
        let dn_st = self.diff_shape_mat_st(s_t_coords);
        let xy = SMatrix::<Dtype, 4, 2>::from([x, y]);
        dn_st * xy
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
    fn calc_k(&self) -> [[Dtype; 8]; 8] {
        println!(
            "\n>>> Calculating Quad2D4N(#{})'s stiffness matrix k{} ......",
            self.id, self.id
        );
        let (ee, nu) = *self.material;
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
        let (ee, nu) = *self.material;
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
    pub fn calc_strain_integration_points(&self) -> [[Dtype; 3]; 4] {
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

    /// Get element node's strain
    pub fn calc_strain_nodes(&self) -> [[Dtype; 3]; 4] {
        let intpt_strain = self.calc_strain_integration_points();
        let mut nodes_strain: [[Dtype; 3]; 4] = [[0.0; 3]; 4];

        let gauss_pt = (1.0 as Dtype) / ((3.0 as Dtype).sqrt());
        let nodes: [[Dtype; 2]; 4] = [[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]];
        let ld_x = -gauss_pt;
        let ld_y = -gauss_pt;
        let rd_x = gauss_pt;
        let rd_y = -gauss_pt;
        let rt_x = gauss_pt;
        let rt_y = gauss_pt;
        let lt_x = -gauss_pt;
        let lt_y = gauss_pt;

        let sigma_xx_intpt = SVector::from([
            intpt_strain[0][0],
            intpt_strain[1][0],
            intpt_strain[2][0],
            intpt_strain[3][0],
        ]);
        let sigma_yy_intpt = SVector::from([
            intpt_strain[0][1],
            intpt_strain[1][1],
            intpt_strain[2][1],
            intpt_strain[3][1],
        ]);
        let sigma_xy_intpt = SVector::from([
            intpt_strain[0][2],
            intpt_strain[1][2],
            intpt_strain[2][2],
            intpt_strain[3][2],
        ]);

        let coord_mat = SMatrix::<Dtype, 4, 4>::from([
            [ld_x * ld_y, ld_x, ld_y, 1.],
            [rd_x * rd_y, rd_x, rd_y, 1.],
            [rt_x * rt_y, rt_x, rt_y, 1.],
            [lt_x * lt_y, lt_x, lt_y, 1.],
        ])
        .transpose();
        let arg_xx: [Dtype; 4] = coord_mat.cholesky().unwrap().solve(&sigma_xx_intpt).into();
        let arg_yy: [Dtype; 4] = coord_mat.cholesky().unwrap().solve(&sigma_yy_intpt).into();
        let arg_xy: [Dtype; 4] = coord_mat.cholesky().unwrap().solve(&sigma_xy_intpt).into();
        print_1darr("arg_xx", &arg_xx, 0., "h");
        print_1darr("arg_yy", &arg_yy, 0., "h");
        print_1darr("arg_xy", &arg_xy, 0., "h");
        nodes_strain
    }

    /// Get element integration points' stress
    pub fn calc_stress_integration_points(&self) -> [[Dtype; 3]; 4] {
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

    /// Output element's calculation result
    pub fn calc_result_info(&self, n_exp: Dtype) -> String {
        let gauss_pt = (1.0 as Dtype) / ((3.0 as Dtype).sqrt());
        let int_pts: [[[Dtype; 2]; 2]; 2] = [
            [[-gauss_pt, -gauss_pt], [gauss_pt, -gauss_pt]],
            [[gauss_pt, gauss_pt], [-gauss_pt, gauss_pt]],
        ];
        format!("\n-----------------------------------------------------------------------------\nElem_Quad2D4N:\n\tId:\t{}\n\tArea: {:-12.6}\n\tMats: {:-12.6} (Young's modulus)\n\t      {:-12.6} (Poisson's ratio)\n\tNodes:{}{}{}{}\n\tStrain: (at integration points)\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\n\tStress: (at integration points)\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\t\t{:-12.6?}\n\n\tStiffness Matrix K{} =  (*10^{})\n{}",
            self.id,
            self.area(),
            self.material.0,
            self.material.1,
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
            self.id(),
            n_exp,
            self.k_string(n_exp),
        )
    }
}

/// Implement zhm::K trait for quadrilateral element
impl<'quad2d4n> K for Quad2D4N<'quad2d4n> {
    /// Cache stiffness matrix for quad element
    fn k(&mut self) -> &CompressedMatrix {
        if self.k_matrix.is_none() {
            self.k_matrix.get_or_insert(compress_matrix(self.calc_k()))
        } else {
            self.k_matrix.as_ref().unwrap()
        }
    }

    /// Print quad element's stiffness matrix
    fn k_printer(&self, n_exp: Dtype) {
        if self.k_matrix.is_none() {
            panic!(
                "!!! Quad2D4N#{}'s k mat is empty! call k() to calc it.",
                self.id
            )
        }

        print!("\nQuad2D4N k{} =  (* 10^{})\n[", self.id, n_exp as i32);
        let elem_stiffness_mat = self.calc_k();
        for row in 0..8 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..8 {
                print!(
                    " {:>-12.6}",
                    elem_stiffness_mat[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                );
            }
            if row == 7 {
                println!("]]");
            } else {
                println!("]");
            }
        }
        print!("\n");
    }

    /// Return quad elem's stiffness matrix's format string
    fn k_string(&self, n_exp: Dtype) -> String {
        if self.k_matrix.is_none() {
            panic!(
                "!!! Quad2D4N#{}'s k mat is empty! call k() to calc it.",
                self.id
            )
        }

        let mut k_matrix = String::new();
        let elem_stiffness_mat = self.calc_k();
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

    /// Get the strain at (xi, eta) inside the element
    fn strain_at_intpt(&mut self) -> Vec<Vec<Dtype>> {
        if self.strain.is_none() {
            self.strain = Some(self.calc_strain_integration_points());
        }
        let intpt_strain = self.strain.unwrap();
        vec![
            intpt_strain[0].to_vec(),
            intpt_strain[1].to_vec(),
            intpt_strain[2].to_vec(),
            intpt_strain[3].to_vec(),
        ]
    }

    /// Get the stress at (xi, eta) inside the element
    fn stress_at_intpt(&mut self) -> Vec<Vec<Dtype>> {
        if self.stress.is_none() {
            self.stress = Some(self.calc_stress_integration_points());
        }
        let intpt_stress = self.stress.unwrap();
        vec![
            intpt_stress[0].to_vec(),
            intpt_stress[1].to_vec(),
            intpt_stress[2].to_vec(),
            intpt_stress[3].to_vec(),
        ]
    }

    /// Get element's info string
    fn info(&self, n_exp: Dtype) -> String {
        self.calc_result_info(n_exp)
    }

    /// Get element's id number
    fn id(&self) -> usize {
        self.id
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
            self.material.0,
            self.material.1,
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
            self.id(),
            self.k_string(0.0)
        )
    }
}
