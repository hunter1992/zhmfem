use crate::{node::Node2D, Dtype, Jacobian2D, K};
use na::*;
use std::fmt::{self, Write};

pub struct Tri2D3N<'tri> {
    pub id: usize,
    pub thick: Dtype,
    pub nodes: [&'tri Node2D; 3],
    pub k_matrix: Option<[[Dtype; 6]; 6]>,
}

impl<'tri> Tri2D3N<'tri> {
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

    /// Calculate the Jacobian matrix of triangle element
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
        let a = [
            xs[idx(0 + 1)] * ys[idx(0 + 2)] - xs[idx(0 + 2)] * ys[idx(0 + 1)],
            xs[idx(1 + 1)] * ys[idx(1 + 2)] - xs[idx(1 + 2)] * ys[idx(1 + 1)],
            xs[idx(2 + 1)] * ys[idx(2 + 2)] - xs[idx(2 + 2)] * ys[idx(2 + 1)],
        ];
        let b = [
            ys[idx(0 + 1)] - ys[idx(0 + 2)],
            ys[idx(1 + 1)] - ys[idx(1 + 2)],
            ys[idx(2 + 1)] - ys[idx(2 + 2)],
        ];
        let c = [
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
        let elasticity_mat = SMatrix::<Dtype, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        let jacobian = Jacobian2D::from(self.jacobian());
        let det_j = jacobian.determinant();

        let b_mat = self.geometry_mat(det_j);
        // Gauss integration, area of standard tri is 0.5
        let core = b_mat.transpose() * elasticity_mat * b_mat * det_j;
        let stiffness_matrix: [[Dtype; 6]; 6] = (0.5 * self.thick * core).into();
        stiffness_matrix
    }

    /// Get element's strain vector
    pub fn calc_strain(&self) -> [Dtype; 3] {
        let jecobian_mat = Jacobian2D::from(self.jacobian());
        let b_mat = self.geometry_mat(jecobian_mat.determinant());
        let elem_nodes_disps = SMatrix::<Dtype, 6, 1>::from(self.disps());
        let strain_vector: [Dtype; 3] = (b_mat * elem_nodes_disps).into();
        strain_vector
    }

    /// Get element's stress vector
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
impl<'tri> K for Tri2D3N<'tri> {
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

pub struct Tri2D6N<'tri> {
    pub id: usize,
    pub thick: Dtype,
    pub nodes: [&'tri Node2D; 3],
    pub k_matrix: Option<[[Dtype; 12]; 12]>,
}

impl<'tri> Tri2D6N<'tri> {
    /// Generate a 2D Tri2D6N element
    pub fn new(id: usize, thick: Dtype, nodes: [&Node2D; 3]) -> Tri2D6N {
        Tri2D6N {
            id,
            thick,
            nodes,
            k_matrix: None,
        }
    }

    /// Get the x-coords of nodes in Tri2D6N element
    pub fn xs(&self) -> [Dtype; 3] {
        let mut x_list = [0.0; 3];
        for i in 0..3 {
            x_list[i] = self.nodes[i].coord[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in Tri2D6N element
    pub fn ys(&self) -> [Dtype; 3] {
        let mut y_list = [0.0; 3];
        for i in 0..3 {
            y_list[i] = self.nodes[i].coord[1];
        }
        y_list
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

    /// Get the 1st row's every single element's cofactor of Determinant A
    ///     | 1  x1  y1 |
    /// A = | 1  x2  y2 |
    ///     | 1  x3  y3 |
    fn a(&self, ith: usize) -> Dtype {
        let xs = self.xs();
        let ys = self.ys();
        let idx = |x: usize| (x % 3) as usize;
        let a = [
            xs[idx(0 + 1)] * ys[idx(0 + 2)] - xs[idx(0 + 2)] * ys[idx(0 + 1)],
            xs[idx(1 + 1)] * ys[idx(1 + 2)] - xs[idx(1 + 2)] * ys[idx(1 + 1)],
            xs[idx(2 + 1)] * ys[idx(2 + 2)] - xs[idx(2 + 2)] * ys[idx(2 + 1)],
        ];
        a[ith]
    }

    /// Get the 2nd row's every single element's cofactor of Determinant A
    ///     | 1  x1  y1 |
    /// A = | 1  x2  y2 |
    ///     | 1  x3  y3 |
    fn b(&self, ith: usize) -> Dtype {
        let ys = self.ys();
        let idx = |x: usize| (x % 3) as usize;
        let b = [
            ys[idx(0 + 1)] - ys[idx(0 + 2)],
            ys[idx(1 + 1)] - ys[idx(1 + 2)],
            ys[idx(2 + 1)] - ys[idx(2 + 2)],
        ];
        b[ith]
    }

    /// Get the 3rd row's every single element's cofactor of Determinant A
    ///     | 1  x1  y1 |
    /// A = | 1  x2  y2 |
    ///     | 1  x3  y3 |
    fn c(&self, ith: usize) -> Dtype {
        let xs = self.xs();
        let idx = |x: usize| (x % 3) as usize;
        let c = [
            xs[idx(0 + 2)] - xs[idx(0 + 1)],
            xs[idx(1 + 2)] - xs[idx(1 + 1)],
            xs[idx(2 + 2)] - xs[idx(2 + 1)],
        ];
        c[ith]
    }

    /// Get area coords Li's value
    /// 自然坐标与面积坐标之间的关系
    /// | L1 |              |a1  b1  c1| | 1 |
    /// | L2 | =  1/(2*A) * |a1  b1  c1|.| x |
    /// | L3 |              |a1  b1  c1| | y |
    fn area_coords(&self, ith: usize) -> impl Fn(Dtype, Dtype) -> Dtype {
        let a = self.a(ith);
        let b = self.b(ith);
        let c = self.c(ith);
        let area: Dtype = self.area();
        move |x: Dtype, y: Dtype| (a + b * x + c * y) / (2.0 * area)
    }

    /// Get shape matrix element N_i
    fn shape_mat_i(&self, ith: usize) -> impl Fn(Dtype, Dtype) -> Dtype {
        /* 形函数矩阵
         * | N1 0  N2 0  N3 0  N4 0  N5 0  N6 0  |
         * | 0  N1 0  N2 0  N3 0  N4 0  N5 0  N6 |
         *
         * N1 = (2L1 - 1) * L1
         * N2 = (2L2 - 1) * L2
         * N3 = (2L3 - 1) * L3
         * N4 = 4 * L1 * L2
         * N5 = 4 * L2 * L3
         * N6 = 4 * L3 * L1
         *
         * 自然坐标与面积坐标之间的关系
         * | L1 |    |a1  b1  c1| | 1  |
         * | L2 | =  |a2  b2  c2|.| x  | * 1/(2*A)
         * | L3 |    |a3  b3  c3| | y  |
         *
         * | 1  |    |1   1   1 | | L1 |
         * | x  | =  |x1  x2  x3|.| L2 |
         * | y  |    |y1  y2  y3| | L2 |
         */
        let n1 = |l1: Dtype, _l4: Dtype| (2.0 * l1 - 1.0) * l1;
        let n2 = |l2: Dtype, _l4: Dtype| (2.0 * l2 - 1.0) * l2;
        let n3 = |l3: Dtype, _l4: Dtype| (2.0 * l3 - 1.0) * l3;
        let n4 = |l1: Dtype, l2: Dtype| 4.0 * l1 * l2;
        let n5 = |l2: Dtype, l3: Dtype| 4.0 * l2 * l3;
        let n6 = |l3: Dtype, l1: Dtype| 4.0 * l3 * l1;

        let shape_mat = [n1, n2, n3, n4, n5, n6];
        shape_mat[ith]
    }

    /// Get B matrix element B_i
    /// | dN/dx    0    |
    /// | 0        dN/dy| =  B
    /// | dN/dy    dN/dx|
    pub fn geometry_mat_i(&self, ith: usize, x_y: (Dtype, Dtype)) -> [[Dtype; 2]; 3] {
        if ith > 5 {
            panic!("!!! Error from Tri2D6N elem func(geometry_mat_i), wrong ith arg.");
        }
        if ith < 3 {
            let l: Dtype = self.area_coords(ith)(x_y.0, x_y.1);
            let arg: Dtype = 0.5 * (4.0 * l - 1.0) / self.area();

            [
                [arg * self.b(ith), 0.0],
                [0.0, arg * self.c(ith)],
                [arg * self.c(ith), arg * self.b(ith)],
            ]
        } else {
            let ith1: usize = ith - 3;
            let ith2: usize = (ith - 2) % 3;
            let l1: Dtype = self.area_coords(ith1)(x_y.0, x_y.1);
            let l2: Dtype = self.area_coords(ith2)(x_y.0, x_y.1);
            let arg1: Dtype = 2.0 * l1 / self.area();
            let arg2: Dtype = 2.0 * l2 / self.area();

            [
                [arg1 * self.b(ith2) + arg2 * self.b(ith1), 0.0],
                [0.0, arg1 * self.c(ith2) + arg2 * self.c(ith1)],
                [
                    arg1 * self.c(ith2) + arg2 * self.c(ith1),
                    arg1 * self.b(ith2) + arg2 * self.b(ith1),
                ],
            ]
        }
    }
}
