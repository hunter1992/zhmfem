use super::triangle::Tri2D3N;
use crate::{compress_matrix, node::Node2D, CompressedMatrix, Data, Dtype, K};
use na::SMatrix;
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
    pub fn get_point_displacement(&self, xi_eta: [Dtype; 3]) -> [Dtype; 2] {
        let s = xi_eta[0];
        let t = xi_eta[1];
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
         * |0    N1   0    N2   0    N3   0    N4|  */
        let x_sign = [-1.0, 1.0, 1.0, -1.0];
        let y_sign = [-1.0, -1.0, 1.0, 1.0];
        move |x: Dtype, y: Dtype| {
            0.25 * (1.0 + x_sign[ith] * x + y_sign[ith] * y + x_sign[ith] * y_sign[ith] * x * y)
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
    fn diff_shape_mat_st(&self, xi_eta: [Dtype; 3]) -> SMatrix<Dtype, 2, 4> {
        let s: Dtype = xi_eta[0];
        let t: Dtype = xi_eta[1];
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
    fn diff_shape_mat_xy(&self, xi_eta: [Dtype; 3]) -> SMatrix<Dtype, 2, 4> {
        let diff_shape_mat_st = self.diff_shape_mat_st(xi_eta);
        let jecobain = self.jacobian(xi_eta);
        jecobain.try_inverse().unwrap() * diff_shape_mat_st
    }

    /// Calculate the Jacobian matrix of quadrilateral
    fn jacobian(&self, xi_eta: [Dtype; 3]) -> SMatrix<Dtype, 2, 2> {
        let x: [Dtype; 4] = self.get_nodes_xcoords();
        let y: [Dtype; 4] = self.get_nodes_ycoords();
        let dn_st = self.diff_shape_mat_st(xi_eta);
        let xy = SMatrix::<Dtype, 4, 2>::from([x, y]);
        dn_st * xy
    }

    /// Geometry matrix B(xi, eta) 用参数坐标表示物理坐标下的几何矩阵
    fn geometry_mat_xy(&self, xi_eta: [Dtype; 3]) -> SMatrix<Dtype, 3, 8> {
        let dn_xy = self.diff_shape_mat_xy(xi_eta);
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
        let gauss_pt = (1.0 as Dtype) / ((3.0 as Dtype).sqrt());
        let int_pts = [
            [[-gauss_pt, -gauss_pt, 0.0], [gauss_pt, -gauss_pt, 0.0]],
            [[gauss_pt, gauss_pt, 0.0], [-gauss_pt, gauss_pt, 0.0]],
        ];

        let mut k_matrix = SMatrix::<Dtype, 8, 8>::from([[0.0; 8]; 8]);
        // 2x2的4点Gauss积分，4个点权重都是1
        for row in 0..2 {
            for col in 0..2 {
                let j = self.jacobian(int_pts[row][col]);
                let det_j = j.determinant();

                let b_mat = self.geometry_mat_xy(int_pts[row][col]);
                let core = b_mat.transpose() * elasticity_mat * b_mat * det_j;
                k_matrix += core;
            }
        }

        //let stiffness_matrix: [[Dtype; 8]; 8] = (self.thick * k_matrix).into();
        let stiffness_matrix: [[Dtype; 8]; 8] = (self.thick * k_matrix).into();
        stiffness_matrix
    }

    /// Get element's strain vector
    fn calc_strain(&self, xi_eta: [Dtype; 3]) -> [Dtype; 3] {
        let b_mat = self.geometry_mat_xy(xi_eta);
        let elem_nodes_disps = SMatrix::<Dtype, 8, 1>::from(self.get_nodes_displacement());
        let strain_vector: [Dtype; 3] = (b_mat * elem_nodes_disps).into();
        strain_vector
    }

    /// Get element's strss vector
    fn calc_stress(&self, xi_eta: [Dtype; 3]) -> [Dtype; 3] {
        let (ee, nu) = *self.material;
        let elasticity_mat = SMatrix::<Dtype, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        let strain = SMatrix::<Dtype, 3, 1>::from(self.calc_strain(xi_eta));
        let stress: [Dtype; 3] = (elasticity_mat * strain).into();
        stress
    }

    /// Get element integration points' strain
    pub fn calc_strain_integration_point(&self) -> [[Dtype; 3]; 4] {
        let mut epsilon: [[Dtype; 3]; 4] = [[0.0; 3]; 4];
        let gauss_pt = (1.0 as Dtype) / ((3.0 as Dtype).sqrt());
        let int_pts: [[Dtype; 3]; 4] = [
            [-gauss_pt, -gauss_pt, 0.0],
            [gauss_pt, -gauss_pt, 0.0],
            [gauss_pt, gauss_pt, 0.0],
            [-gauss_pt, gauss_pt, 0.0],
        ];
        for idx in 0..4 {
            epsilon[idx] = self.calc_strain(int_pts[idx]);
        }
        epsilon
    }

    /// Get element integration points' stress
    pub fn calc_stress_integration_point(&self) -> [[Dtype; 3]; 4] {
        let mut sigma: [[Dtype; 3]; 4] = [[0.0; 3]; 4];
        let gauss_pt = (1.0 as Dtype) / ((3.0 as Dtype).sqrt());
        let int_pts: [[Dtype; 3]; 4] = [
            [-gauss_pt, -gauss_pt, 0.0],
            [gauss_pt, -gauss_pt, 0.0],
            [gauss_pt, gauss_pt, 0.0],
            [-gauss_pt, gauss_pt, 0.0],
        ];
        for idx in 0..4 {
            sigma[idx] = self.calc_stress(int_pts[idx]);
        }
        sigma
    }

    /// Print element's strain value
    pub fn print_strain(&self, xi_eta: [Dtype; 3]) {
        let strain = self.calc_strain(xi_eta);
        println!(
            "\nelem[{}] strain:\n\tE_xx = {:-16.6}\n\tE_yy = {:-16.6}\n\tE_xy = {:-16.6}",
            self.id, strain[0], strain[1], strain[2]
        );
    }

    /// Print element's stress value
    pub fn print_stress(&self, xi_eta: [Dtype; 3]) {
        let stress = self.calc_stress(xi_eta);
        println!(
            "\nelem[{}] stress:\n\tS_xx = {:-16.6}\n\tS_yy = {:-16.6}\n\tS_xy = {:-16.6}",
            self.id, stress[0], stress[1], stress[2]
        );
    }

    /// Output element's calculation result
    pub fn calc_result_info(&self, n_exp: Dtype) -> String {
        let gauss_pt = (1.0 as Dtype) / ((3.0 as Dtype).sqrt());
        let int_pts: [[[Dtype; 3]; 2]; 2] = [
            [[-gauss_pt, -gauss_pt, 0.0], [gauss_pt, -gauss_pt, 0.0]],
            [[gauss_pt, gauss_pt, 0.0], [-gauss_pt, gauss_pt, 0.0]],
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
    fn strain_intpt(&mut self) -> Data {
        if self.strain.is_none() {
            self.strain
                .get_or_insert(self.calc_strain_integration_point());
        }
        let mut data: Vec<[Dtype; 3]> = Vec::with_capacity(4);
        for idx in 0..4 {
            data.push(self.strain.unwrap()[idx]);
        }
        Data::Dim2(data)
    }

    /// Get the stress at (xi, eta) inside the element
    fn stress_intpt(&mut self) -> Data {
        if self.stress.is_none() {
            self.stress
                .get_or_insert(self.calc_stress_integration_point());
        }
        let mut data: Vec<[Dtype; 3]> = Vec::with_capacity(4);
        for idx in 0..4 {
            data.push(self.stress.unwrap()[idx]);
        }
        Data::Dim2(data)
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
        let int_pts: [[[Dtype; 3]; 2]; 2] = [
            [[-gauss_pt, -gauss_pt, 0.0], [gauss_pt, -gauss_pt, 0.0]],
            [[gauss_pt, gauss_pt, 0.0], [-gauss_pt, gauss_pt, 0.0]],
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
