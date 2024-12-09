use crate::{node::Node2D, Dtype, Jacobian2D, K};
use na::SMatrix;
use std::fmt::{self, Write};

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

    /// Get shape matrix element N_i
    fn shape_mat_i(&self, ith: usize) -> impl Fn(Dtype, Dtype) -> Dtype {
        /* The shape mat of tri elem:
         * |N1  0   N2  0   N3  0 |
         * |0   N1  0   N2  0   N3|
         * 输入参数i用于选择输出哪个Ni
         */
        let area = self.area();
        let xs = self.get_nodes_xcoords();
        let ys = self.get_nodes_ycoords();
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

    /// Get any point's disps vector in tri element
    pub fn point_disp(&self, point_coord: [Dtype; 2]) -> [Dtype; 2] {
        let x = point_coord[0];
        let y = point_coord[1];
        let n0 = self.shape_mat_i(0usize)(x, y);
        let n1 = self.shape_mat_i(1usize)(x, y);
        let n2 = self.shape_mat_i(2usize)(x, y);

        let disps = self.get_nodes_displacement();
        let u = n0 * disps[0] + n1 * disps[2] + n2 * disps[4];
        let v = n0 * disps[1] + n1 * disps[3] + n2 * disps[5];
        [u, v]
    }

    /// Calculate the Jacobian matrix of tri2d3n element
    fn jacobian(&self) -> [[Dtype; 2]; 2] {
        let x: [Dtype; 3] = self.get_nodes_xcoords();
        let y: [Dtype; 3] = self.get_nodes_ycoords();
        let dx_21 = x[1] - x[0];
        let dx_31 = x[2] - x[0];
        let dy_21 = y[1] - y[0];
        let dy_31 = y[2] - y[0];
        [[dx_21, dx_31], [dy_21, dy_31]]
    }

    /// Element's B matrix, B mat is the combination of diff(N)
    fn geometry_mat(&self, det_j: Dtype) -> SMatrix<Dtype, 3, 6> {
        let x: [Dtype; 3] = self.get_nodes_xcoords();
        let y: [Dtype; 3] = self.get_nodes_ycoords();
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

        let b_mat = self.geometry_mat(det_j);
        // Gauss integration, area of standard tri is 0.5
        let core = b_mat.transpose() * elasticity_mat * b_mat * det_j;
        let stiffness_matrix: [[Dtype; 6]; 6] = (0.5 * self.thick * core).into();
        stiffness_matrix
    }

    /// Get element's strain vector, the strain in CST elem is a const
    fn calc_strain(&self) -> [Dtype; 3] {
        let jecobian_mat = Jacobian2D::from(self.jacobian()).transpose();
        let b_mat = self.geometry_mat(jecobian_mat.determinant());
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
