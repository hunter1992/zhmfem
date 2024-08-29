//use crate::{matrix_block_fill, matrix_hstack, print_2darr};
use crate::{node::Node2D, Dtype, Jacobian2D, K};
use na::*;
use std::fmt::{self, Write};

pub struct Tri2D3N<'tri2d3n> {
    pub id: usize,
    pub thick: Dtype,
    pub nodes: [&'tri2d3n Node2D; 3],
    pub k_matrix: Option<[[Dtype; 6]; 6]>,
}

impl<'tri2d3n> Tri2D3N<'tri2d3n> {
    /// Generate a 2D Tri2D3N element
    pub fn new(id: usize, thick: Dtype, nodes: [&Node2D; 3]) -> Tri2D3N {
        Tri2D3N {
            id,
            thick,
            nodes,
            k_matrix: None,
        }
    }

    /// Get the x-coords of nodes in tri element
    pub fn xs(&self) -> [Dtype; 3] {
        let mut x_list = [0.0; 3];
        for i in 0..3 {
            x_list[i] = self.nodes[i].coord[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in tri element
    pub fn ys(&self) -> [Dtype; 3] {
        let mut y_list = [0.0; 3];
        for i in 0..3 {
            y_list[i] = self.nodes[i].coord[1];
        }
        y_list
    }

    /// Get triangle element area value
    pub fn area(&self) -> Dtype {
        let x = self.xs();
        let y = self.ys();
        let dx_21 = x[1] - x[0];
        let dx_31 = x[2] - x[0];
        let dy_21 = y[1] - y[0];
        let dy_31 = y[2] - y[0];
        0.5 * (dx_21 * dy_31 - dx_31 * dy_21).abs()
    }

    /// Calculate the Jacobian matrix of tri2d3n element
    pub fn jacobian(&self) -> [[Dtype; 2]; 2] {
        let x: [Dtype; 3] = self.xs();
        let y: [Dtype; 3] = self.ys();
        let dx_21 = x[1] - x[0];
        let dx_31 = x[2] - x[0];
        let dy_21 = y[1] - y[0];
        let dy_31 = y[2] - y[0];
        [[dx_21, dx_31], [dy_21, dy_31]]
    }

    /// Get nodes' disps vector in tri element
    pub fn disps(&self) -> [Dtype; 6] {
        let mut disps = [0.0; 6];
        for idx in 0..3 {
            disps[2 * idx] = *self.nodes[idx].disps[0].borrow();
            disps[2 * idx + 1] = *self.nodes[idx].disps[1].borrow();
        }
        disps
    }

    /// Get nodes's force vector in tri element
    pub fn forces(&self) -> [Dtype; 6] {
        let mut forces = [0.0; 6];
        for idx in 0..3 {
            forces[2 * idx] = *self.nodes[idx].forces[0].borrow();
            forces[2 * idx + 1] = *self.nodes[idx].forces[1].borrow();
        }
        forces
    }

    /// Get any point's disps vector in tri element
    pub fn point_disp(&self, point_coord: [Dtype; 2]) -> [Dtype; 2] {
        let x = point_coord[0];
        let y = point_coord[1];
        let n0 = self.shape_mat_i(0usize)(x, y);
        let n1 = self.shape_mat_i(1usize)(x, y);
        let n2 = self.shape_mat_i(2usize)(x, y);

        let disps = self.disps();
        let u = n0 * disps[0] + n1 * disps[2] + n2 * disps[4];
        let v = n0 * disps[1] + n1 * disps[3] + n2 * disps[5];
        [u, v]
    }

    /// Get shape matrix element N_i
    fn shape_mat_i(&self, ith: usize) -> impl Fn(Dtype, Dtype) -> Dtype {
        /* The shape mat of tri elem:
         * |N1  0   N2  0   N3  0 |
         * |0   N1  0   N2  0   N3|
         * 输入参数i用于选择输出哪个Ni
         */
        let area = self.area();
        let xs = self.xs();
        let ys = self.ys();
        let idx = |x: usize| (x % 3) as usize;
        let a: [Dtype; 3] = [
            xs[idx(0 + 1)] * ys[idx(0 + 2)] - xs[idx(0 + 2)] * ys[idx(0 + 1)],
            xs[idx(1 + 1)] * ys[idx(1 + 2)] - xs[idx(1 + 2)] * ys[idx(1 + 1)],
            xs[idx(2 + 1)] * ys[idx(2 + 2)] - xs[idx(2 + 2)] * ys[idx(2 + 1)],
        ];
        let b: [Dtype; 3] = [
            ys[idx(0 + 1)] - ys[idx(0 + 2)],
            ys[idx(1 + 1)] - ys[idx(1 + 2)],
            ys[idx(2 + 1)] - ys[idx(2 + 2)],
        ];
        let c: [Dtype; 3] = [
            xs[idx(0 + 2)] - xs[idx(0 + 1)],
            xs[idx(1 + 2)] - xs[idx(1 + 1)],
            xs[idx(2 + 2)] - xs[idx(2 + 1)],
        ];
        move |x: Dtype, y: Dtype| 0.5 * (a[ith] + x * b[ith] + y * c[ith]) / area
    }

    /// Element's B matrix, B mat is the combination of diff(N)
    fn geometry_mat(&self, det_j: Dtype) -> SMatrix<Dtype, 3, 6> {
        let x: [Dtype; 3] = self.xs();
        let y: [Dtype; 3] = self.ys();
        let h_mat = SMatrix::<Dtype, 3, 4>::from([
            [y[2] - y[0], 0.0, x[0] - x[2]],
            [y[0] - y[1], 0.0, x[1] - x[0]],
            [0.0, x[0] - x[2], y[2] - y[0]],
            [0.0, x[1] - x[0], y[0] - y[1]],
        ]) / (det_j.abs());

        let q_mat = SMatrix::<Dtype, 4, 6>::from([
            [-1., -1., 0., 0.],
            [0., 0., -1., -1.],
            [1., 0., 0., 0.],
            [0., 0., 1., 0.],
            [0., 1., 0., 0.],
            [0., 0., 0., 1.],
        ]);

        h_mat * q_mat
    }

    /// Calculate element stiffness matrix K
    /// Return a 6x6 matrix, elements are Dtype
    fn calc_k(&self, material_args: (Dtype, Dtype)) -> [[Dtype; 6]; 6] {
        println!(
            "\n>>> Calculating Tri2D3N(#{})'s stiffness matrix k{} ......",
            self.id, self.id
        );
        let (ee, nu) = material_args;
        let elasticity_mat = (ee / (1.0 - nu * nu))
            * (SMatrix::<Dtype, 3, 3>::from([
                [1.0, nu, 0.0],
                [nu, 1.0, 0.0],
                [0.0, 0.0, 0.5 * (1.0 - nu)],
            ])
            .transpose());

        let jacobian = Jacobian2D::from(self.jacobian()).transpose();
        let det_j = jacobian.determinant();

        let b_mat = self.geometry_mat(det_j);
        // Gauss integration, area of standard tri is 0.5
        let core = b_mat.transpose() * elasticity_mat * b_mat * det_j;
        let stiffness_matrix: [[Dtype; 6]; 6] = (0.5 * self.thick * core).into();
        stiffness_matrix
    }

    /// Get element's strain vector, the strain in CST elem is a const
    pub fn calc_strain(&self) -> [Dtype; 3] {
        let jecobian_mat = Jacobian2D::from(self.jacobian()).transpose();
        let b_mat = self.geometry_mat(jecobian_mat.determinant());
        let elem_nodes_disps = SMatrix::<Dtype, 6, 1>::from(self.disps());
        let strain_vector: [Dtype; 3] = (b_mat * elem_nodes_disps).into();
        strain_vector
    }

    /// Get element's stress vector, the stress in CST elem is a const
    pub fn calc_stress(&self, material_args: (Dtype, Dtype)) -> [Dtype; 3] {
        let (ee, nu) = material_args;
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
            "\nelem[{}] strain:\n\tE_xx = {:-9.6}\n\tE_yy = {:-9.6}\n\tE_xy = {:-9.6}",
            self.id, strain[0], strain[1], strain[2]
        );
    }

    /// Print element's stress value
    pub fn print_stress(&self, material_args: (Dtype, Dtype)) {
        let stress = self.calc_stress(material_args);
        println!(
            "\nelem[{}] stress:\n\tS_xx = {:-9.6}\n\tS_yy = {:-9.6}\n\tS_xy = {:-9.6}",
            self.id, stress[0], stress[1], stress[2]
        );
    }
}

/// Implement zhm::K trait for triangle element
impl<'tri2d3n> K for Tri2D3N<'tri2d3n> {
    type Kmatrix = [[Dtype; 6]; 6];

    /// Cache stiffness matrix for triangle element
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
                    " {:>-10.6} ",
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
                    " {:>-10.6} ",
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
    fn stress(&self, _xyz: [Dtype; 3], material: (Dtype, Dtype)) -> Vec<Dtype> {
        self.calc_stress(material).to_vec()
    }

    /// Get element's info string
    fn info(&self) -> String {
        format!("\n--------------------------------------------------------------------\nElement_2D Info:\n\tId:     {}\n\tArea:   {}\n\tType:   Tri2D3N
\tNodes: {}\n\t       {}\n\t       {}\n",
            self.id,
            self.area(),
            self.nodes[0],
            self.nodes[1],
            self.nodes[2]
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
            "\nElement_2D Info:\n\tId:     {}\n\tArea:   {}\n\tType:   Tri2D3N
\tNodes: {}\n\t       {}\n\t       {}",
            self.id,
            self.area(),
            self.nodes[0],
            self.nodes[1],
            self.nodes[2]
        )
    }
}

pub struct Tri2D6N<'tri2d6n> {
    pub id: usize,
    pub thick: Dtype,
    pub nodes: [&'tri2d6n Node2D; 6],
    pub k_matrix: Option<[[Dtype; 12]; 12]>,
}

impl<'tri2d6n> Tri2D6N<'tri2d6n> {
    /// Generate a 2D Tri2D6N element
    pub fn new(id: usize, thick: Dtype, nodes: [&Node2D; 6]) -> Tri2D6N {
        Tri2D6N {
            id,
            thick,
            nodes,
            k_matrix: None,
        }
    }

    /// Get the x-coords of nodes in Tri2D6N element
    pub fn xs(&self) -> [Dtype; 6] {
        let mut x_list = [0.0; 6];
        for i in 0..6 {
            x_list[i] = self.nodes[i].coord[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in Tri2D6N element
    pub fn ys(&self) -> [Dtype; 6] {
        let mut y_list = [0.0; 6];
        for i in 0..6 {
            y_list[i] = self.nodes[i].coord[1];
        }
        y_list
    }

    /// Get nodes' disps vector in tri2d6n element
    pub fn disps(&self) -> [Dtype; 12] {
        let mut disps = [0.0; 12];
        for idx in 0..6 {
            disps[2 * idx] = *self.nodes[idx].disps[0].borrow();
            disps[2 * idx + 1] = *self.nodes[idx].disps[1].borrow();
        }
        disps
    }

    /// Get nodes' force vector in tri2d6n element
    pub fn forces(&self) -> [Dtype; 12] {
        let mut forces = [0.0; 12];
        for idx in 0..6 {
            forces[2 * idx] = *self.nodes[idx].disps[0].borrow();
            forces[2 * idx + 1] = *self.nodes[idx].disps[1].borrow();
        }
        forces
    }

    /// Get Tri2D6N element area value
    pub fn area(&self) -> Dtype {
        let x = self.xs();
        let y = self.ys();
        let dx_21 = x[1] - x[0];
        let dx_31 = x[2] - x[0];
        let dy_21 = y[1] - y[0];
        let dy_31 = y[2] - y[0];
        0.5 * (dx_21 * dy_31 - dx_31 * dy_21).abs()
    }

    /// Get shape matrix element N_i(s, t)
    /// shape matrix:
    /// | N1  0   N2  0   N3  0   N4  0   N5  0   N6  0  |
    /// | 0   N1  0   N2  0   N3  0   N4  0   N5  0   N6 |
    ///
    /// N1 = (1 - s - t) * (1 - 2s - 2t)
    /// N2 = s * (2s - 1)
    /// N3 = t * (2t - 1)
    /// N4 = 4 * s * (1 - s - t)
    /// N5 = 4 * s * t
    /// N6 = 4 * t * (1 - s - t)
    /// for every s and t in [0, 1]
    #[inline]
    pub fn shape_mat(&self, s: Dtype, t: Dtype) -> [[Dtype; 12]; 2] {
        let n_mat: [[Dtype; 12]; 2] = [
            [
                (1.0 - s - t) * (1.0 - 2.0 * s - 2.0 * t),
                0.0,
                s * (2.0 * s - 1.0),
                0.0,
                t * (2.0 * t - 1.0),
                0.0,
                4.0 * s * (1.0 - s - t),
                0.0,
                4.0 * s * t,
                0.0,
                4.0 * t * (1.0 - s - t),
                0.0,
            ],
            [
                0.0,
                (1.0 - s - t) * (1.0 - 2.0 * s - 2.0 * t),
                0.0,
                s * (2.0 * s - 1.0),
                0.0,
                t * (2.0 * t - 1.0),
                0.0,
                4.0 * s * (1.0 - s - t),
                0.0,
                4.0 * s * t,
                0.0,
                4.0 * t * (1.0 - s - t),
            ],
        ];
        n_mat
    }

    /// Return B mat of Tri2D6N element
    pub fn geometry_mat_xy(&self, s_t: &[Dtype]) -> [[Dtype; 12]; 3] {
        let s = s_t[0];
        let t = s_t[1];
        let j = self.jacobian(s_t);
        let j_det = j[1][1] * j[0][0] - j[0][1] * j[1][0];
        let j_inv = Matrix2::from([[j[1][1], -j[0][1]], [j[1][0], j[0][0]]]).transpose();

        let dn1_ds: Dtype = 4.0 * (s + t) - 3.0;
        let dn2_ds: Dtype = 4.0 * s - 1.0;
        let dn3_ds: Dtype = 0.0;
        let dn4_ds: Dtype = 4.0 * (1.0 - t) - 8.0 * s;
        let dn5_ds: Dtype = 4.0 * t;
        let dn6_ds: Dtype = -4.0 * t;

        let dn1_dt: Dtype = 4.0 * (s + t) - 3.0;
        let dn2_dt: Dtype = 0.0;
        let dn3_dt: Dtype = 4.0 * t - 1.0;
        let dn4_dt: Dtype = 4.0 * (s + t) - 3.0;
        let dn5_dt: Dtype = 4.0 * (s + t) - 3.0;
        let dn6_dt: Dtype = 4.0 * (s + t) - 3.0;

        let dn_ds = (SMatrix::<Dtype, 6, 2>::from([
            [dn1_ds, dn2_ds, dn3_ds, dn4_ds, dn5_ds, dn6_ds],
            [dn1_dt, dn2_dt, dn3_dt, dn4_dt, dn5_dt, dn6_dt],
        ])
        .transpose())
            / j_det;
        let dn_dx: [[Dtype; 6]; 2] = ((j_inv * dn_ds) / j_det).transpose().into();
        let b_mat: [[Dtype; 12]; 3] = [
            [
                dn_dx[0][0],
                0.0,
                dn_dx[0][1],
                0.0,
                dn_dx[0][2],
                0.0,
                dn_dx[0][3],
                0.0,
                dn_dx[0][4],
                0.0,
                dn_dx[0][5],
                0.0,
            ],
            [
                dn_dx[1][0],
                0.0,
                dn_dx[1][1],
                0.0,
                dn_dx[1][2],
                0.0,
                dn_dx[1][3],
                0.0,
                dn_dx[1][4],
                0.0,
                dn_dx[1][5],
                0.0,
            ],
            [
                dn_dx[1][0],
                dn_dx[0][0],
                dn_dx[1][1],
                dn_dx[0][1],
                dn_dx[1][2],
                dn_dx[0][2],
                dn_dx[1][3],
                dn_dx[0][3],
                dn_dx[1][4],
                dn_dx[0][4],
                dn_dx[1][5],
                dn_dx[0][5],
            ],
        ];
        b_mat
    }

    /// Calculate the Jacobian matrix of tri2d6n element
    fn jacobian(&self, s_t: &[Dtype]) -> [[Dtype; 2]; 2] {
        let s = s_t[0];
        let t = s_t[1];
        let x: [Dtype; 6] = self.xs();
        let y: [Dtype; 6] = self.ys();

        let dn1_ds = |s: Dtype, t: Dtype| 4.0 * (s + t) - 3.0;
        let dn2_ds = |s: Dtype, _t: Dtype| 4.0 * s - 1.0;
        let dn3_ds = |_s: Dtype, _t: Dtype| 0.0;
        let dn4_ds = |s: Dtype, t: Dtype| 4.0 * (1.0 - t) - 8.0 * s;
        let dn5_ds = |_s: Dtype, t: Dtype| 4.0 * t;
        let dn6_ds = |_s: Dtype, t: Dtype| -4.0 * t;
        let dn_ds = [dn1_ds, dn2_ds, dn3_ds, dn4_ds, dn5_ds, dn6_ds];

        let dn1_dt = |s: Dtype, t: Dtype| 4.0 * (s + t) - 3.0;
        let dn2_dt = |_s: Dtype, _t: Dtype| 0.0;
        let dn3_dt = |_s: Dtype, t: Dtype| 4.0 * t - 1.0;
        let dn4_dt = |s: Dtype, t: Dtype| 4.0 * (s + t) - 3.0;
        let dn5_dt = |s: Dtype, t: Dtype| 4.0 * (s + t) - 3.0;
        let dn6_dt = |s: Dtype, t: Dtype| 4.0 * (s + t) - 3.0;
        let dn_dt = [dn1_dt, dn2_dt, dn3_dt, dn4_dt, dn5_dt, dn6_dt];

        let mut dx_ds: Dtype = 0.0;
        let mut dx_dt: Dtype = 0.0;
        let mut dy_ds: Dtype = 0.0;
        let mut dy_dt: Dtype = 0.0;
        for i in 0..6 {
            dx_ds += x[i] * dn_ds[i](s, t);
            dx_dt += x[i] * dn_dt[i](s, t);
            dy_ds += y[i] * dn_ds[i](s, t);
            dy_dt += y[i] * dn_dt[i](s, t);
        }
        let j_mat = [[dx_ds, dx_dt], [dy_ds, dy_dt]];

        if (dx_ds * dy_dt - dx_dt * dy_ds) < 0.0 {
            panic!(">>> !!! The determinant of Jacobian matrix is negative !!! \n>>>     from src/elem/dim2/linear/fn geometry_mat_xy");
        } else {
            j_mat
        }
    }

    /// Calculate element stiffness matrix K
    /// Return a 12x12 matrix, elements are Dtype
    pub fn calc_k(&self, material_args: (Dtype, Dtype)) -> [[Dtype; 12]; 12] {
        println!(
            "\n>>> Calculating Tri2D3N(#{})'s stiffness matrix k{} ......",
            self.id, self.id
        );
        let thick = self.thick;
        let (ee, nu) = material_args;
        let elasticity_mat = SMatrix::<Dtype, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        let det_j = |s_t: &[Dtype]| {
            Jacobian2D::from(self.jacobian(s_t))
                .transpose()
                .determinant()
        };
        let int_point_2 = [
            [0.6666666667, 0.1666666667],
            [0.1666666667, 0.6666666667],
            [0.1666666667, 0.1666666667],
        ];
        let weight = [0.3333333333, 0.3333333333, 0.3333333333];

        let b_mat_1 =
            SMatrix::<Dtype, 12, 3>::from(self.geometry_mat_xy(&int_point_2[0])).transpose();
        let b_mat_2 =
            SMatrix::<Dtype, 12, 3>::from(self.geometry_mat_xy(&int_point_2[1])).transpose();
        let b_mat_3 =
            SMatrix::<Dtype, 12, 3>::from(self.geometry_mat_xy(&int_point_2[2])).transpose();
        let core = b_mat_1.transpose()
            * elasticity_mat
            * b_mat_1
            * det_j(&int_point_2[0])
            * weight[0]
            + b_mat_2.transpose() * elasticity_mat * b_mat_2 * det_j(&int_point_2[1]) * weight[1]
            + b_mat_3.transpose() * elasticity_mat * b_mat_3 * det_j(&int_point_2[2]) * weight[2];

        let stiffness_matrix: [[Dtype; 12]; 12] = (0.5 * thick * core).into();
        stiffness_matrix
    }

    /// Get element's strain vector, the strain in CST elem is a const
    pub fn calc_strain(&self, xyz: [Dtype; 3]) -> [Dtype; 3] {
        todo!()
    }

    /// Get element's stress vector, the stress in CST elem is a const
    pub fn calc_stress(&self, xyz: [Dtype; 3], material_args: (Dtype, Dtype)) -> [Dtype; 3] {
        todo!()
    }

    /// Print element's strain value
    pub fn print_strain(&self, xyz: [Dtype; 3]) {
        let strain = self.calc_strain(xyz);
        println!(
            "\nelem[{}], strain @ ({:-6.3}, {:-6.3}):\n\tE_xx = {:-9.6}\n\tE_yy = {:-9.6}\n\tE_xy = {:-9.6}",
            self.id, &xyz[0], &xyz[1], strain[0], strain[1], strain[2]
        );
    }

    /// Print element's stress value
    pub fn print_stress(&self, xyz: [Dtype; 3], material_args: (Dtype, Dtype)) {
        let stress = self.calc_stress(xyz, material_args);
        println!(
            "\nelem[{}], stress @ ({:-6.3}, {:-6.3}):\n\tS_xx = {:-9.6}\n\tS_yy = {:-9.6}\n\tS_xy = {:-9.6}",
            self.id, &xyz[0], &xyz[1], stress[0], stress[1], stress[2]
        );
    }
}

/// Implement zhm::K trait for 2D 6nodes triangle
impl<'tri2d6n> K for Tri2D6N<'tri2d6n> {
    type Kmatrix = [[Dtype; 12]; 12];

    /// Cache stiffness matrix for triangle element
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

    /// Print triangle element's stiffness matrix
    fn k_printer(&self, n_exp: Dtype) {
        if self.k_matrix.is_none() {
            panic!(
                "!!! Tri2D3N#{}'s k mat is empty! call k() to calc it.",
                self.id
            );
        }

        print!("\nTri2D3N k{} =  (* 10^{})\n[", self.id, n_exp as u8);
        for row in 0..12 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..12 {
                print!(
                    " {:>-10.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                );
            }
            if row == 11 {
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
        for row in 0..12 {
            if row == 0 {
                write!(k_matrix, "[[").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, " [").expect("!!! Write tri k_mat failed!");
            }
            for col in 0..12 {
                write!(
                    k_matrix,
                    " {:>-10.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                )
                .expect("!!! Write tri k_mat failed!");
            }
            if row == 11 {
                write!(k_matrix, "]]").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, "]\n").expect("!!! Write tri k_mat failed!");
            }
        }
        k_matrix
    }

    /// Get the strain at (x,y) inside the element
    /// In linear triangle elem, strain is a const
    fn strain(&self, xyz: [Dtype; 3]) -> Vec<Dtype> {
        self.calc_strain(xyz).to_vec()
    }

    /// Get the stress at (x,y) inside the element
    /// In linear triangle elem, stress is a const
    fn stress(&self, xyz: [Dtype; 3], material: (Dtype, Dtype)) -> Vec<Dtype> {
        self.calc_stress(xyz, material).to_vec()
    }

    /// Get element's info string
    fn info(&self) -> String {
        format!("\n--------------------------------------------------------------------\nElement_2D Info:\n\tId:     {}\n\tArea:   {}\n\tType:   Tri2D6N
\tNodes: {}\n\t       {}\n\t       {}\n\t       {}\n\t       {}\n\t       {}\n\t",
            self.id,
            self.area(),
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
            self.nodes[3],
            self.nodes[4],
            self.nodes[5]
        )
    }

    /// Get element's id number
    fn id(&self) -> usize {
        self.id
    }
}

impl fmt::Display for Tri2D6N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nElement_2D Info:\n\tId:     {}\n\tArea:   {}\n\tType:   Tri2D3N
\tNodes: {}\n\t       {}\n\t       {}\n\t       {}\n\t       {}\n\t       {}\n\t",
            self.id,
            self.area(),
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
            self.nodes[3],
            self.nodes[4],
            self.nodes[5]
        )
    }
}
