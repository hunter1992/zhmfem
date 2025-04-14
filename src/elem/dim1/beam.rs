use crate::{node::Node2D, Dtype, K, SS};
use na::SMatrix;
use std::fmt::{self, Write};

/// 细长梁以及纯弯梁，抽象成具有如下特征的理论模型：
///     1) 仅用x坐标刻画所有位置(即细长之含义)；
///     2) 变形主要为垂直于x轴的挠度
///     3) 节点受载荷为横力和弯矩
/// 但是，节点的自由度不为1，包含：挠度(y)、转角(theta)，
/// 此亦是梁单元为非协调单元的原因。
pub struct Beam1D2N<'beam1d2n> {
    pub id: usize,
    pub moment_of_inertia: Dtype,
    pub cross_sectional_area: Dtype,
    pub nodes: [&'beam1d2n Node2D; 2],
    pub strain: Option<SS>,
    pub stress: Option<SS>,
    pub k_matrix: Option<[[Dtype; 4]; 4]>,
    pub material: &'beam1d2n (Dtype, Dtype),
}

impl<'beam1d2n> Beam1D2N<'beam1d2n> {
    /// Generate a Beam2D2N element
    pub fn new(
        id: usize,
        moment_of_inertia: Dtype,
        cross_sectional_area: Dtype,
        nodes: [&'beam1d2n Node2D; 2],
        material: &'beam1d2n (Dtype, Dtype),
    ) -> Self {
        Beam1D2N {
            id,
            moment_of_inertia,
            cross_sectional_area,
            nodes,
            strain: None,
            stress: None,
            k_matrix: None,
            material,
        }
    }

    /// Set element material_args
    pub fn set_material(&mut self, material_args: &'beam1d2n (Dtype, Dtype)) {
        self.material = material_args;
    }

    /// Get difference in x-coordinates of rod nodes
    pub fn dx(&self) -> Dtype {
        let x = self.get_nodes_xcoords();
        x[1] - x[0]
    }

    /// Get difference in y-coordinates of rod nodes
    pub fn dy(&self) -> Dtype {
        let y = self.get_nodes_ycoords();
        y[1] - y[0]
    }

    /// Get Beam2D2N element length
    pub fn length(&self) -> Dtype {
        let dx = self.dx();
        let dy = self.dy();
        (dx * dx + dy * dy).sqrt()
    }

    /// Get the x-coords of nodes in tri element
    pub fn get_nodes_xcoords(&self) -> [Dtype; 2] {
        let mut x_list = [0.0; 2];
        for i in 0..2 {
            x_list[i] = self.nodes[i].coords[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in tri element
    pub fn get_nodes_ycoords(&self) -> [Dtype; 2] {
        let mut y_list = [0.0; 2];
        for i in 0..2 {
            y_list[i] = self.nodes[i].coords[1];
        }
        y_list
    }

    /// Get nodes' disps vector in beam2d2n element
    pub fn get_nodes_displacement(&self) -> [Dtype; 4] {
        let mut disps = [0.0; 4];
        for idx in 0..2 {
            disps[2 * idx] = self.nodes[idx].displs.borrow()[0];
            disps[2 * idx + 1] = self.nodes[idx].displs.borrow()[1];
        }
        disps
    }

    /// Get nodes's force vector in tri element
    pub fn get_nodes_force(&self) -> [Dtype; 4] {
        let mut forces = [0.0; 4];
        for idx in 0..2 {
            forces[2 * idx] = self.nodes[idx].forces.borrow()[0];
            forces[2 * idx + 1] = self.nodes[idx].forces.borrow()[1];
        }
        forces
    }

    /// Get shape matrix element N_i
    /// The shape mat of Beam1D2N element: [N1 N2 N3 N4]
    /// which is a 1x4 mat.    \xi = x/L
    /// N1 = 1*(1 + 0*\xi - 3*\xi^2 + 2*\xi^3)
    /// N2 = L*(0 + 1*\xi - 2*\xi^2 + 1*\xi^3)
    /// N3 = 1*(0 + 0*\xi + 3*\xi^2 - 2*\xi^3)
    /// N4 = L*(0 + 0*\xi - 1*\xi^2 + 1*\xi^3)
    fn shape_func_st(&self, ith: usize) -> impl Fn(Dtype) -> Dtype {
        let l1 = self.length();
        let l2 = l1 * l1;
        let l3 = l1 * l1 * l1;
        let col1: [Dtype; 4] = [1.0, 0.0, 0.0, 0.0];
        let col2: [Dtype; 4] = [0.0, l1, 0.0, 0.0];
        let col3: [Dtype; 4] = [-3.0, -2.0 * l1, 3.0, -l1];
        let col4: [Dtype; 4] = [2.0, l1, -2.0, l1];
        move |x: Dtype| {
            col1[ith] + col2[ith] * x / l1 + col3[ith] * x * x / l2 + col4[ith] * x * x * x / l3
        }
    }

    /// Get any point's disps vector in beam element
    /// Return a 1x2 matrix: [deflection, rotation]
    pub fn point_disp(&self, point_coord: [Dtype; 1]) -> [Dtype; 2] {
        let x = point_coord[0];
        let n0 = self.shape_func_st(0usize)(x);
        let n1 = self.shape_func_st(1usize)(x);
        let n2 = self.shape_func_st(2usize)(x);
        let n3 = self.shape_func_st(3usize)(x);

        let disps = self.get_nodes_displacement();
        let v = n0 * disps[0] + n2 * disps[2];
        let theta = n1 * disps[1] + n3 * disps[3];
        [v, theta]
    }

    /// Element's B matrix, B mat is the combination of diff(N)
    fn geometry_mat(&self, xi: Dtype) -> SMatrix<Dtype, 1, 4> {
        let l: Dtype = self.length();
        SMatrix::<Dtype, 1, 4>::from([
            [12.0 * xi - 6.0],
            [(6.0 * xi - 4.0) * l],
            [-12.0 * xi + 6.0],
            [(6.0 * xi - 2.0) * l],
        ])
    }

    /// Return a 4x4 matrix, elements are Dtype
    fn calc_k(&self) -> [[Dtype; 4]; 4] {
        println!(
            "\n>>> Calculating Beam2D2N(#{})'s local stiffness matrix k{} ......",
            self.id, self.id
        );
        let moi = self.moment_of_inertia;
        let ee = self.material.0;
        let l = self.length();
        let l3 = l * l * l;
        let gauss_pt: Dtype = (3.0_f32).sqrt() / 6.0;
        let int_pts: [Dtype; 2] = [0.5 - gauss_pt, 0.5 + gauss_pt];

        let mut k_matrix = SMatrix::<Dtype, 4, 4>::from([[0.0; 4]; 4]);
        for idx in 0..2 {
            let b_mat = self.geometry_mat(int_pts[idx]);
            let core = 0.5 * b_mat.transpose() * b_mat;
            k_matrix = k_matrix + core;
        }

        let stiffness_matrix: [[Dtype; 4]; 4] = (ee * moi * k_matrix / l3).into();
        stiffness_matrix
    }

    /// Get the strain at point (x) inside the element
    fn calc_strain(&self, xyz: [Dtype; 3]) -> [Dtype; 1] {
        let x = xyz[0];
        let y = xyz[1];
        let xi = x / self.length();
        let b_mat = self.geometry_mat(xi);
        let elem_nodes_disps = SMatrix::<Dtype, 4, 1>::from(self.get_nodes_displacement());
        let strain: [Dtype; 1] = (-y * b_mat * elem_nodes_disps).into();
        strain
    }

    /// Get element's stress vector, the stress in CST elem is a const
    fn calc_stress(&self, xyz: [Dtype; 3]) -> [Dtype; 1] {
        let ee = self.material.0;
        [ee * self.calc_strain(xyz)[0]]
    }

    /// Get the strain at integration point
    fn calc_strain_integration_point(&self) -> [[Dtype; 1]; 2] {
        let mut epsilon: [[Dtype; 1]; 2] = [[0.0; 1]; 2];
        let gauss_pt: Dtype = (3.0_f32).sqrt() / 6.0;
        let int_pts: [Dtype; 2] = [0.5 - gauss_pt, 0.5 + gauss_pt];
        for idx in 0..2 {
            epsilon[idx] = self.calc_strain([int_pts[idx], 0.0, 0.0]);
        }
        epsilon
    }

    /// Get the stress at integration point
    fn calc_stress_integration_point(&self) -> [[Dtype; 1]; 2] {
        let mut sigma: [[Dtype; 1]; 2] = [[0.0; 1]; 2];
        let gauss_pt: Dtype = (3.0_f32).sqrt() / 6.0;
        let int_pts: [Dtype; 2] = [0.5 - gauss_pt, 0.5 + gauss_pt];
        for idx in 0..2 {
            sigma[idx] = self.calc_strain([int_pts[idx], 0.0, 0.0]);
        }
        sigma
    }

    /// Print element's strain value
    pub fn print_strain(&self, xyz: [Dtype; 3]) {
        let strain = self.calc_strain(xyz);
        println!(
            "\nelem[{}] strain:\n\tE_xx = {:-16.6}\n\tE_yy = {:-16.6}\n\tE_xy = {:-16.6}",
            self.id, strain[0], 0.0, 0.0
        );
    }

    /// Print element's stress value
    pub fn print_stress(&self, xyz: [Dtype; 3]) {
        let stress = self.calc_stress(xyz);
        println!(
            "\nelem[{}] stress:\n\tS_xx = {:-16.6}\n\tS_yy = {:-16.6}\n\tS_xy = {:-16.6}",
            self.id, stress[0], 0.0, 0.0
        );
    }
}

/// Implement zhm::K trait for triangle element
impl<'beam2d2n> K for Beam1D2N<'beam2d2n> {
    type Kmatrix = [[Dtype; 4]; 4];

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
                "!!! Beam2D2N#{}'s k mat is empty! call k() to calc it.",
                self.id
            );
        }

        print!("\nBeam2D2N k{} =  (* 10^{})\n[", self.id, n_exp as u8);
        for row in 0..4 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..4 {
                print!(
                    " {:>-12.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                );
            }
            if row == 3 {
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
        for row in 0..4 {
            if row == 0 {
                write!(k_matrix, "[[").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, " [").expect("!!! Write tri k_mat failed!");
            }
            for col in 0..4 {
                write!(
                    k_matrix,
                    " {:>-12.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                )
                .expect("!!! Write tri k_mat failed!");
            }
            if row == 3 {
                write!(k_matrix, "]]").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, "]\n").expect("!!! Write tri k_mat failed!");
            }
        }
        k_matrix
    }

    /// Get the strain at (x,y) inside the element
    fn strain_intpt(&mut self) -> &SS {
        if self.strain.is_none() {
            self.strain
                .get_or_insert(SS::Dim1(self.calc_strain([0.5, 0.0, 0.0])))
        } else {
            self.strain.as_ref().unwrap()
        }
    }

    /// Get the stress at (x,y) inside the element
    fn stress_intpt(&mut self) -> &SS {
        if self.stress.is_none() {
            self.stress
                .get_or_insert(SS::Dim1(self.calc_stress([0.5, 0.0, 0.0])))
        } else {
            self.stress.as_ref().unwrap()
        }
    }

    /// Get element's info string
    fn info(&self, n_exp: Dtype) -> String {
        format!("\n-----------------------------------------------------------------------------\nElem_Beam2D2N:\n\tId:\t{}\n\tArea: {:>-12.4}\n\tMoI : {:>-12.4} (Moment of inertia)\n\tMats: {:>-12.4} (Young's modulus)\n\t      {:>-12.4} (Poisson's ratio)\n\tNodes:{}{}\n\tStiffness Matrix K{} = \n{}",
            self.id,
            self.cross_sectional_area,
            self.moment_of_inertia,
            self.material.0,
            self.material.1,
            self.nodes[0],
            self.nodes[1],
            self.id(),
            self.k_string(n_exp),
        )
    }

    /// Get element's id number
    fn id(&self) -> usize {
        self.id
    }
}

impl fmt::Display for Beam1D2N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nElem_Tri2D3N:\n\tId:\t{}\n\tArea: {:>-12.4}\n\tMats: {:>-12.4} (Young's modulus)\n\t      {:>-12.4} (Poisson's ratio)\n\tNodes:{}{}{}",
            self.id,
            self.cross_sectional_area,
            self.moment_of_inertia,
            self.material.0,
            self.material.1,
            self.nodes[0],
            self.nodes[1],
        )
    }
}
