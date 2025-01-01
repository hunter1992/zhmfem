use crate::{node::Node2D, Dtype, Jacobian2D, K};
use na::SMatrix;
use std::fmt::{self, Write};

pub struct Tri2D6N<'tri2d6n> {
    pub id: usize,
    pub thick: Dtype,
    pub nodes: [&'tri2d6n Node2D; 6],
    pub k_matrix: Option<[[Dtype; 12]; 12]>,
    pub material: &'tri2d6n (Dtype, Dtype),
}

impl<'tri2d6n> Tri2D6N<'tri2d6n> {
    /// Generate a 2D Tri2D6N element
    pub fn new(
        id: usize,
        thick: Dtype,
        nodes: [&'tri2d6n Node2D; 6],
        material: &'tri2d6n (Dtype, Dtype),
    ) -> Self {
        Tri2D6N {
            id,
            thick,
            nodes,
            k_matrix: None,
            material,
        }
    }

    /// Set element material_args
    pub fn set_material(&mut self, material_args: &'tri2d6n (Dtype, Dtype)) {
        self.material = material_args;
    }

    /// Get triangle element area value
    pub fn area(&self) -> Dtype {
        let x = self.get_nodes_xcoords();
        let y = self.get_nodes_ycoords();
        let dx_21 = x[1] - x[0];
        let dx_31 = x[2] - x[0];
        let dy_21 = y[1] - y[0];
        let dy_31 = y[2] - y[0];
        0.5 * (dx_21 * dy_31 - dx_31 * dy_21).abs()
    }

    /// Get the x-coords of nodes in tri element
    pub fn get_nodes_xcoords(&self) -> [Dtype; 6] {
        let mut x_list = [0.0; 6];
        for i in 0..6 {
            x_list[i] = self.nodes[i].coords[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in tri element
    pub fn get_nodes_ycoords(&self) -> [Dtype; 6] {
        let mut y_list = [0.0; 6];
        for i in 0..6 {
            y_list[i] = self.nodes[i].coords[1];
        }
        y_list
    }

    /// Get nodes' disps vector in tri element
    pub fn get_nodes_displacement(&self) -> [Dtype; 12] {
        let mut disps = [0.0; 12];
        for idx in 0..6 {
            disps[2 * idx] = self.nodes[idx].displs.borrow()[0];
            disps[2 * idx + 1] = self.nodes[idx].displs.borrow()[1];
        }
        disps
    }

    /// Get nodes's force vector in tri element
    pub fn get_nodes_force(&self) -> [Dtype; 12] {
        let mut forces = [0.0; 12];
        for idx in 0..6 {
            forces[2 * idx] = self.nodes[idx].forces.borrow()[0];
            forces[2 * idx + 1] = self.nodes[idx].forces.borrow()[1];
        }
        forces
    }

    /// a,b,c通过计算下式的代数余子式获得：
    /// | 1  x1  y1 |     | a1  b1  c1 |
    /// | 1  x2  y2 |     | a2  b2  c2 |
    /// | 1  x3  y3 |     | a3  b3  c3 |
    /// 对应关系是右边ai或bi或ci位置上的代数余子式
    fn coef_a(&self, ith: usize) -> Dtype {
        let xs = self.get_nodes_xcoords();
        let ys = self.get_nodes_ycoords();
        let idx = |x: usize| (x % 3) as usize;
        let a: [Dtype; 3] = [
            xs[idx(0 + 1)] * ys[idx(0 + 2)] - xs[idx(0 + 2)] * ys[idx(0 + 1)],
            xs[idx(1 + 1)] * ys[idx(1 + 2)] - xs[idx(1 + 2)] * ys[idx(1 + 1)],
            xs[idx(2 + 1)] * ys[idx(2 + 2)] - xs[idx(2 + 2)] * ys[idx(2 + 1)],
        ];
        a[ith]
    }

    /// Calculate coefficient b
    fn coef_b(&self, ith: usize) -> Dtype {
        let ys = self.get_nodes_ycoords();
        let idx = |x: usize| (x % 3) as usize;
        let b: [Dtype; 3] = [
            ys[idx(0 + 1)] - ys[idx(0 + 2)],
            ys[idx(1 + 1)] - ys[idx(1 + 2)],
            ys[idx(2 + 1)] - ys[idx(2 + 2)],
        ];
        b[ith]
    }

    /// Calculate coefficient c
    fn coef_c(&self, ith: usize) -> Dtype {
        let xs = self.get_nodes_xcoords();
        let idx = |x: usize| (x % 3) as usize;
        let c: [Dtype; 3] = [
            xs[idx(0 + 2)] - xs[idx(0 + 1)],
            xs[idx(1 + 2)] - xs[idx(1 + 1)],
            xs[idx(2 + 2)] - xs[idx(2 + 1)],
        ];
        c[ith]
    }

    /// Transform physical coordinates into parametric coordinates
    /// | L1 |                  | a1  b1  c1 |   | 1 |
    /// | L2 |  =  (1 / 2*area) | a2  b2  c2 | * | x |
    /// | L3 |                  | a3  b3  c3 |   | y |
    fn coords_transform_xtol(&self, xyz: [Dtype; 3]) -> [Dtype; 3] {
        let x = xyz[0];
        let y = xyz[1];
        let area = self.area();
        let mut l: [Dtype; 3] = [0.0; 3];
        for idx in 0..3 {
            l[idx] = 0.5 * (self.coef_a(idx) + self.coef_b(idx) * x + self.coef_c(idx) * y) / area;
        }
        l
    }

    /// Transform parametric coordinates into physical coordinates
    /// | 1 |     |  1   1   1 |   | L1 |
    /// | x |  =  | x1  x2  x3 | * | L2 |
    /// | y |     | y1  y2  y3 |   | L3 |
    fn coords_transform_ltox(&self, l123: [Dtype; 3]) -> [Dtype; 3] {
        let l1 = l123[0];
        let l2 = l123[1];
        let l3 = l123[2];
        let mut xyz: [Dtype; 3] = [0.0; 3];
        let x1 = self.get_nodes_xcoords()[0];
        let x2 = self.get_nodes_xcoords()[1];
        let x3 = self.get_nodes_xcoords()[2];
        let y1 = self.get_nodes_ycoords()[0];
        let y2 = self.get_nodes_ycoords()[1];
        let y3 = self.get_nodes_ycoords()[2];
        xyz[0] = x1 * l1 + x2 * l2 + x3 * l3;
        xyz[1] = y1 * l1 + y2 * l2 + y3 * l3;
        xyz
    }

    /// Get shape matrix element N_i(s, t)
    /// The shape mat of tri2d6n elem:
    /// |N1  0   N2  0   N3  0   N4  0   N5  0   N6  0 |
    /// |0   N1  0   N2  0   N3  0   N4  0   N5  0   N6|
    ///
    /// N1 = (2L1 - 1) * L1
    /// N2 = (2L2 - 1) * L2
    /// N3 = (2L3 - 1) * L3
    /// N4 = 4 * L1 * L2
    /// N5 = 4 * L2 * L3
    /// N6 = 4 * L3 * L1
    ///
    /// 输入参数i用于选择输出哪个Ni
    fn shape_mat_i(&self, ith: usize) -> impl Fn(Dtype, Dtype, Dtype) -> Dtype {
        let n1 = move |l1: Dtype, _: Dtype, _: Dtype| (2.0 * l1 - 1.0) * l1;
        let n2 = move |_: Dtype, l2: Dtype, _: Dtype| (2.0 * l2 - 1.0) * l2;
        let n3 = move |_: Dtype, _: Dtype, l3: Dtype| (2.0 * l3 - 1.0) * l3;
        let n4 = move |l1: Dtype, l2: Dtype, _: Dtype| 4.0 * l1 * l2;
        let n5 = move |_: Dtype, l2: Dtype, l3: Dtype| 4.0 * l2 * l3;
        let n6 = move |l1: Dtype, _: Dtype, l3: Dtype| 4.0 * l3 * l1;
        let n = [n1, n2, n3, n4, n5, n6];
        n[ith]
    }

    /// point displacement inner tri2d6n
    pub fn point_disp(&self, point_coord_local: [Dtype; 3]) -> [Dtype; 2] {
        let [l1, l2, l3] = self.coords_transform_xtol(point_coord_local);
        let n0 = self.shape_mat_i(0usize)(l1, l2, l3);
        let n1 = self.shape_mat_i(1usize)(l1, l2, l3);
        let n2 = self.shape_mat_i(2usize)(l1, l2, l3);
        let n3 = self.shape_mat_i(3usize)(l1, l2, l3);
        let n4 = self.shape_mat_i(4usize)(l1, l2, l3);
        let n5 = self.shape_mat_i(5usize)(l1, l2, l3);

        let disps = self.get_nodes_displacement();
        let u = n0 * disps[0]
            + n1 * disps[2]
            + n2 * disps[4]
            + n3 * disps[6]
            + n4 * disps[8]
            + n5 * disps[10];
        let v = n0 * disps[1]
            + n1 * disps[3]
            + n2 * disps[5]
            + n3 * disps[7]
            + n4 * disps[9]
            + n5 * disps[11];
        [u, v]
    }
}

pub struct Tri2D3N<'tri2d3n> {
    pub id: usize,
    pub thick: Dtype,
    pub nodes: [&'tri2d3n Node2D; 3],
    pub k_matrix: Option<[[Dtype; 6]; 6]>,
    pub material: &'tri2d3n (Dtype, Dtype),
}

impl<'tri2d3n> Tri2D3N<'tri2d3n> {
    /// Generate a 2D Tri2D3N element
    pub fn new(
        id: usize,
        thick: Dtype,
        nodes: [&'tri2d3n Node2D; 3],
        material: &'tri2d3n (Dtype, Dtype),
    ) -> Self {
        Tri2D3N {
            id,
            thick,
            nodes,
            k_matrix: None,
            material,
        }
    }

    /// Set element material_args
    pub fn set_material(&mut self, material_args: &'tri2d3n (Dtype, Dtype)) {
        self.material = material_args;
    }

    /// Get triangle element area value
    pub fn area(&self) -> Dtype {
        let x = self.get_nodes_xcoords();
        let y = self.get_nodes_ycoords();
        let dx_21 = x[1] - x[0];
        let dx_31 = x[2] - x[0];
        let dy_21 = y[1] - y[0];
        let dy_31 = y[2] - y[0];
        0.5 * (dx_21 * dy_31 - dx_31 * dy_21).abs()
    }

    /// Get the x-coords of nodes in tri element
    pub fn get_nodes_xcoords(&self) -> [Dtype; 3] {
        let mut x_list = [0.0; 3];
        for i in 0..3 {
            x_list[i] = self.nodes[i].coords[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in tri element
    pub fn get_nodes_ycoords(&self) -> [Dtype; 3] {
        let mut y_list = [0.0; 3];
        for i in 0..3 {
            y_list[i] = self.nodes[i].coords[1];
        }
        y_list
    }

    /// Get nodes' disps vector in tri element
    pub fn get_nodes_displacement(&self) -> [Dtype; 6] {
        let mut disps = [0.0; 6];
        for idx in 0..3 {
            disps[2 * idx] = self.nodes[idx].displs.borrow()[0];
            disps[2 * idx + 1] = self.nodes[idx].displs.borrow()[1];
        }
        disps
    }

    /// Get nodes's force vector in tri element
    pub fn get_nodes_force(&self) -> [Dtype; 6] {
        let mut forces = [0.0; 6];
        for idx in 0..3 {
            forces[2 * idx] = self.nodes[idx].forces.borrow()[0];
            forces[2 * idx + 1] = self.nodes[idx].forces.borrow()[1];
        }
        forces
    }

    /// Coefficient a of shape function
    /// u(x, y) = A0 + A1 * x + A2 * y
    /// v(x, y) = B0 + B1 * x + B2 * y
    ///
    /// A0 = (a1*u1 + a2*u2 + a3*u3) / 2A
    /// A1 = (b1*u1 + b2*u2 + b3*u3) / 2A
    /// A2 = (c1*u1 + c2*u2 + c3*u3) / 2A
    ///
    /// B0 = (a1*v1 + a2*v2 + a3*v3) / 2A
    /// B1 = (b1*v1 + b2*v2 + b3*v3) / 2A
    /// B2 = (c1*v1 + c2*v2 + c3*v3) / 2A
    ///
    /// a,b,c通过计算下式的代数余子式获得：
    /// | 1  x1  y1 |     | a1  b1  c1 |
    /// | 1  x2  y2 |     | a2  b2  c2 |
    /// | 1  x3  y3 |     | a3  b3  c3 |
    /// 对应关系是右边ai或bi或ci位置上的代数余子式
    fn coef_a(&self, ith: usize) -> Dtype {
        let xs = self.get_nodes_xcoords();
        let ys = self.get_nodes_ycoords();
        let idx = |x: usize| (x % 3) as usize;
        let a: [Dtype; 3] = [
            xs[idx(0 + 1)] * ys[idx(0 + 2)] - xs[idx(0 + 2)] * ys[idx(0 + 1)],
            xs[idx(1 + 1)] * ys[idx(1 + 2)] - xs[idx(1 + 2)] * ys[idx(1 + 1)],
            xs[idx(2 + 1)] * ys[idx(2 + 2)] - xs[idx(2 + 2)] * ys[idx(2 + 1)],
        ];
        a[ith]
    }

    /// calculate coefficient b
    fn coef_b(&self, ith: usize) -> Dtype {
        let ys = self.get_nodes_ycoords();
        let idx = |x: usize| (x % 3) as usize;
        let b: [Dtype; 3] = [
            ys[idx(0 + 1)] - ys[idx(0 + 2)],
            ys[idx(1 + 1)] - ys[idx(1 + 2)],
            ys[idx(2 + 1)] - ys[idx(2 + 2)],
        ];
        b[ith]
    }

    /// calculate coefficient c
    fn coef_c(&self, ith: usize) -> Dtype {
        let xs = self.get_nodes_xcoords();
        let idx = |x: usize| (x % 3) as usize;
        let c: [Dtype; 3] = [
            xs[idx(0 + 2)] - xs[idx(0 + 1)],
            xs[idx(1 + 2)] - xs[idx(1 + 1)],
            xs[idx(2 + 2)] - xs[idx(2 + 1)],
        ];
        c[ith]
    }

    /// Get Tri2D3N's shape matrix element N_i
    /// 线性三角形单元的三个插值函数就是三个面积坐标
    /// The shape mat of Tri2D3N:
    /// |N1  0   N2  0   N3  0 |
    /// |0   N1  0   N2  0   N3|
    /// 即Ni = Li (i=1,2,3)
    /// 输入参数i用于选择输出哪个Ni (i=1,2,3)
    fn shape_mat_i(&self, ith: usize) -> impl Fn(Dtype, Dtype) -> Dtype {
        let area = self.area();
        let a = self.coef_a(ith);
        let b = self.coef_b(ith);
        let c = self.coef_c(ith);
        move |x: Dtype, y: Dtype| 0.5 * (a + x * b + y * c) / area
    }

    /// Get any point's disps vector in tri element
    pub fn point_disp(&self, point_coord_local: [Dtype; 2]) -> [Dtype; 2] {
        let s = point_coord_local[0];
        let t = point_coord_local[1];
        let n0 = self.shape_mat_i(0usize)(s, t);
        let n1 = self.shape_mat_i(1usize)(s, t);
        let n2 = self.shape_mat_i(2usize)(s, t);

        let disps = self.get_nodes_displacement();
        let u = n0 * disps[0] + n1 * disps[2] + n2 * disps[4];
        let v = n0 * disps[1] + n1 * disps[3] + n2 * disps[5];
        [u, v]
    }

    /// Calculate the Jacobian matrix of tri2d3n element
    /// J = [[dx/ds, dy/ds]
    ///      [dx/dt, dy/dt]]
    fn jacobian(&self) -> [[Dtype; 2]; 2] {
        let x: [Dtype; 3] = self.get_nodes_xcoords();
        let y: [Dtype; 3] = self.get_nodes_ycoords();
        let dx_21 = x[1] - x[0];
        let dx_31 = x[2] - x[0];
        let dy_21 = y[1] - y[0];
        let dy_31 = y[2] - y[0];
        [[dx_21, dx_31], [dy_21, dy_31]]
    }

    /// Calculate H matrix
    /// The H matrix transforms the strain in the parameter coordinate system to
    /// the physical coordinate system
    ///
    /// epsilon = Du/Dx
    ///         = Du/Ds * Ds/Dx
    ///         = J^(-1) * Du/Ds
    ///         = H * Epsilon
    ///
    /// |     du/dx     |        | du/ds |
    /// |     dv/dy     |  = H * | du/dt |
    /// | du/dy + dv/dx |        | dv/ds |
    ///                          | dv/dt |
    ///
    /// epsilon: Strain in physical coordinate system
    /// Epsilon: Strain in parametric coordinate system
    fn geometry_transform_mat(&self, j_mat: [[Dtype; 2]; 2]) -> SMatrix<Dtype, 3, 4> {
        let det_j: Dtype = j_mat[0][0] * j_mat[1][1] - j_mat[0][1] * j_mat[1][0];
        if det_j < 0.0 {
            panic!("!!! Det of Jacobian < 0 !");
        }
        let h = SMatrix::<Dtype, 3, 4>::from([
            [j_mat[1][1], 0.0, -j_mat[1][0]],
            [-j_mat[0][1], 0.0, j_mat[0][0]],
            [0.0, -j_mat[1][0], j_mat[1][1]],
            [0.0, j_mat[0][0], -j_mat[0][1]],
        ]);
        h / det_j
    }

    /// Calculate Q matrix
    /// Q matrix is the  differentiation of N(s, t):
    ///
    /// | du/ds |   | d /ds,   0   |   | u |
    /// | du/dt | = | d /dt,   0   | * |   |
    /// | dv/ds |   |   0  , d /ds |   | v |
    /// | dv/dt |   |   0  , d /dt |   |   |
    ///                                                             | u1 |
    ///             | d /ds,   0   |                                | v1 |
    ///           = | d /dt,   0   | * | N1, 0 , N2, 0 , N3, 0  | * | u2 |
    ///             |   0  , d /ds |   | 0 , N1, 0 , N2, 0 , N3 |   | v2 |
    ///             |   0  , d /dt |                                | u3 |
    ///                                                             | v3 |
    ///           = Q * | u1 |
    ///                 | v1 |
    ///                 | u2 |
    ///                 | v2 |
    ///                 | u3 |
    ///                 | v3 |
    ///
    /// Q = | d /ds,   0   |
    ///     | d /dt,   0   | * | N1, 0 , N2, 0 , N3, 0  |
    ///     |   0  , d /ds |   | 0 , N1, 0 , N2, 0 , N3 |
    ///     |   0  , d /dt |
    ///
    ///   = | b1  0  b2  0  b3  0 |
    ///     | c1  0  c2  0  c3  0 |
    ///     |  0  b2  0  b3  0 b1 |
    ///     |  0  b2  0  b3  0 b1 |
    fn diff_shape_mat(&self) -> SMatrix<Dtype, 4, 6> {
        let dn1_ds: Dtype = self.coef_b(0);
        let dn1_dt: Dtype = self.coef_c(0);
        let dn2_ds: Dtype = self.coef_b(1);
        let dn2_dt: Dtype = self.coef_c(1);
        let dn3_ds: Dtype = self.coef_b(2);
        let dn3_dt: Dtype = self.coef_c(2);
        SMatrix::<Dtype, 4, 6>::from([
            [dn1_ds, dn1_dt, 0.0, 0.0],
            [0.0, 0.0, dn1_ds, dn1_dt],
            [dn2_ds, dn2_dt, 0.0, 0.0],
            [0.0, 0.0, dn2_ds, dn2_dt],
            [dn3_ds, dn3_dt, 0.0, 0.0],
            [0.0, 0.0, dn3_ds, dn3_dt],
        ])
    }

    /// Element's B matrix is the combination of diff(N)
    fn geometry_mat(&self) -> SMatrix<Dtype, 3, 6> {
        let h_mat = self.geometry_transform_mat(self.jacobian());
        let q_mat = self.diff_shape_mat();
        h_mat * q_mat
    }

    /// Return a 6x6 matrix, elements are Dtype
    fn calc_k(&self) -> [[Dtype; 6]; 6] {
        println!(
            "\n>>> Calculating Tri2D3N(#{})'s local stiffness matrix k{} ......",
            self.id, self.id
        );
        let (ee, nu) = *self.material;
        let elasticity_mat = (ee / (1.0 - nu * nu))
            * (SMatrix::<Dtype, 3, 3>::from([
                [1.0, nu, 0.0],
                [nu, 1.0, 0.0],
                [0.0, 0.0, 0.5 * (1.0 - nu)],
            ])
            .transpose());

        let jacobian = Jacobian2D::from(self.jacobian()).transpose();
        let det_j = jacobian.determinant();
        let b_mat = self.geometry_mat();

        // Gauss integration, area of standard tri is 0.5
        let core = b_mat.transpose() * elasticity_mat * b_mat * det_j;
        let stiffness_matrix: [[Dtype; 6]; 6] = (0.5 * self.thick * core).into();
        stiffness_matrix
    }

    /// Get element's strain vector, the strain in CST elem is a const
    fn calc_strain(&self) -> [Dtype; 3] {
        let b_mat = self.geometry_mat();
        let elem_nodes_disps = SMatrix::<Dtype, 6, 1>::from(self.get_nodes_displacement());
        let strain_vector: [Dtype; 3] = (b_mat * elem_nodes_disps).into();
        strain_vector
    }

    /// Get element's stress vector, the stress in CST elem is a const
    fn calc_stress(&self) -> [Dtype; 3] {
        let (ee, nu) = *self.material;
        let elasticity_mat = SMatrix::<Dtype, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        let strain = SMatrix::<Dtype, 3, 1>::from(self.calc_strain());
        let stress: [Dtype; 3] = (elasticity_mat * strain).into();
        stress
    }

    /// Print element's strain value
    pub fn print_strain(&self) {
        let strain = self.calc_strain();
        println!(
            "\nelem[{}] strain:\n\tE_xx = {:-16.6}\n\tE_yy = {:-16.6}\n\tE_xy = {:-16.6}",
            self.id, strain[0], strain[1], strain[2]
        );
    }

    /// Print element's stress value
    pub fn print_stress(&self) {
        let stress = self.calc_stress();
        println!(
            "\nelem[{}] stress:\n\tS_xx = {:-16.6}\n\tS_yy = {:-16.6}\n\tS_xy = {:-16.6}",
            self.id, stress[0], stress[1], stress[2]
        );
    }
}

/// Implement zhm::K trait for triangle element
impl<'tri2d3n> K for Tri2D3N<'tri2d3n> {
    type Kmatrix = [[Dtype; 6]; 6];

    /// Cache stiffness matrix for triangle element
    fn k(&mut self) -> &Self::Kmatrix
    where
        Self::Kmatrix: std::ops::Index<usize>,
    {
        if self.k_matrix.is_none() {
            self.k_matrix.get_or_insert(self.calc_k())
        } else {
            self.k_matrix.as_ref().unwrap()
        }
    }

    /// Print triangle element's stiffness matrix
    fn k_printer(&self, n_exp: Dtype) {
        if self.k_matrix.is_none() {
            panic!(
                "!!! Tri2D3N#{}'s k mat is empty! call k() to calc it.",
                self.id
            );
        }

        print!("\nTri2D3N k{} =  (* 10^{})\n[", self.id, n_exp as u8);
        for row in 0..6 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..6 {
                print!(
                    " {:>-12.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                );
            }
            if row == 5 {
                println!("]]");
            } else {
                println!("]");
            }
        }
        print!("\n");
    }

    /// Return triangle elem's stiffness matrix's format string
    fn k_string(&self, n_exp: Dtype) -> String {
        let mut k_matrix = String::new();
        for row in 0..6 {
            if row == 0 {
                write!(k_matrix, "[[").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, " [").expect("!!! Write tri k_mat failed!");
            }
            for col in 0..6 {
                write!(
                    k_matrix,
                    " {:>-12.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                )
                .expect("!!! Write tri k_mat failed!");
            }
            if row == 5 {
                write!(k_matrix, "]]").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, "]\n").expect("!!! Write tri k_mat failed!");
            }
        }
        k_matrix
    }

    /// Get the strain at (x,y) inside the element
    /// In linear triangle elem, strain is a const
    fn strain(&self, _xyz: [Dtype; 3]) -> Vec<Dtype> {
        self.calc_strain().to_vec()
    }

    /// Get the stress at (x,y) inside the element
    /// In linear triangle elem, stress is a const
    fn stress(&self, _xyz: [Dtype; 3]) -> Vec<Dtype> {
        self.calc_stress().to_vec()
    }

    /// Get element's info string
    fn info(&self, n_exp: Dtype) -> String {
        format!("\n-----------------------------------------------------------------------------\nElem_Tri2D3N:\n\tId:\t{}\n\tArea: {:>-12.4}\n\tMats: {:>-12.4} (Young's modulus)\n\t      {:>-12.4} (Poisson's ratio)\n\tNodes:{}{}{}\n\tStrain:\n\t\t{:-12.6?}\n\n\tStress:\n\t\t{:-12.6?}\n\n\tStiffness Matrix K{} = \n{}",
            self.id,
            self.area(),
            self.material.0,
            self.material.1,
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
            self.strain([0.0, 0.0, 0.0]),
            self.stress([0.0, 0.0, 0.0]),
            self.id(),
            self.k_string(n_exp),
        )
    }

    /// Get element's id number
    fn id(&self) -> usize {
        self.id
    }
}

impl fmt::Display for Tri2D3N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nElem_Tri2D3N:\n\tId:\t{}\n\tArea: {:>-12.4}\n\tMats: {:>-12.4} (Young's modulus)\n\t      {:>-12.4} (Poisson's ratio)\n\tNodes:{}{}{}",
            self.id,
            self.area(),
            self.material.0,
            self.material.1,
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
        )
    }
}
