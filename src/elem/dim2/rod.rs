use crate::data::{CompressedMatrix, Data, Dtype};
use crate::node::Node2D;
use crate::port::K;
use crate::tool::{compress_matrix, print_2darr};
use na::SMatrix;
use std::fmt::{self, Write};

pub struct Rod2D2N<'rod2d2n> {
    pub id: usize,
    pub cross_sectional_area: Dtype,
    pub nodes: [&'rod2d2n Node2D; 2],
    pub strain: Option<[Dtype; 1]>,
    pub stress: Option<[Dtype; 1]>,
    pub k_matrix: Option<CompressedMatrix>,
    pub material: &'rod2d2n (Dtype, Dtype),
}

impl<'rod2d2n> Rod2D2N<'rod2d2n> {
    /// Generate a 2D Rod2D2N element
    pub fn new(
        id: usize,
        cross_sectional_area: Dtype,
        nodes: [&'rod2d2n Node2D; 2],
        material: &'rod2d2n (Dtype, Dtype),
    ) -> Self {
        Rod2D2N {
            id,
            cross_sectional_area,
            nodes,
            strain: None,
            stress: None,
            k_matrix: None,
            material,
        }
    }

    /// Set element material_args
    pub fn set_material(&mut self, material_args: &'rod2d2n (Dtype, Dtype)) {
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

    /// Get Rod2D2N element length
    pub fn length(&self) -> Dtype {
        let dx = self.dx();
        let dy = self.dy();
        (dx * dx + dy * dy).sqrt()
    }

    /// Cosine of the angle between the rod and the horizontal axis
    pub fn cosine(&self) -> Dtype {
        ((self.dx() as f64) / self.length() as f64) as Dtype
    }

    /// Sine of the angle between the rod and the horizontal axis
    pub fn sine(&self) -> Dtype {
        ((self.dy() as f64) / self.length() as f64) as Dtype
    }

    /// Get the x-coords of nodes in Rod1D2N element
    pub fn get_nodes_xcoords(&self) -> [Dtype; 2] {
        let mut x_list = [0.0; 2];
        for i in 0..2 {
            x_list[i] = self.nodes[i].coords[0];
        }
        x_list
    }

    /// Get the x-coords of nodes in Rod1D2N element
    pub fn get_nodes_ycoords(&self) -> [Dtype; 2] {
        let mut y_list = [0.0; 2];
        for i in 0..2 {
            y_list[i] = self.nodes[i].coords[1];
        }
        y_list
    }

    /// Get transformation matrix, trans global to local
    pub fn trans_mat(&self) -> SMatrix<Dtype, 2, 4> {
        let cos = self.cosine();
        let sin = self.sine();
        let mat = [[cos, 0.0], [sin, 0.0], [0.0, cos], [0.0, sin]];
        SMatrix::<Dtype, 2, 4>::from(mat)
    }

    /// Get nodes' disps vector in Rod1D2N element
    pub fn get_nodes_displacement(&self) -> [Dtype; 4] {
        let mut disps = [0.0; 4];
        for idx in 0..2 {
            disps[idx * 2] = self.nodes[idx].displs.borrow()[0];
            disps[idx * 2 + 1] = self.nodes[idx].displs.borrow()[1];
        }
        disps
    }

    /// Get nodes's force vector in Rod1D2N element
    pub fn get_nodes_force(&self) -> [Dtype; 4] {
        let mut forces = [0.0; 4];
        for idx in 0..2 {
            forces[idx * 2] = self.nodes[idx].forces.borrow()[0];
            forces[idx * 2 + 1] = self.nodes[idx].forces.borrow()[1];
        }
        forces
    }

    /// Get axial force in Rod2D2N element
    pub fn axial_force(&self) -> Dtype {
        self.calc_stress()[0] * self.cross_sectional_area
    }

    /// Get shape matrix element N_i
    /// The shape mat of rod elem:
    /// [N1 N2]  which is a 1x2 mat
    /// N1 = 1 - x/L = 1 - epsilon
    /// N2 = x/L     =     epsilon
    fn shape_mat_i(&self, ith: usize) -> impl Fn(Dtype) -> Dtype {
        let a: [Dtype; 2] = [1.0, 0.0];
        let b: [Dtype; 2] = [-1.0, 1.0];
        let length = self.length();
        move |x: Dtype| a[ith] + b[ith] * x / length
    }

    /// Get point's disp in element coord
    pub fn point_disp(&self, point_coord: [Dtype; 1]) -> [Dtype; 1] {
        let x = point_coord[0];
        let n0 = self.shape_mat_i(0usize)(x);
        let n1 = self.shape_mat_i(1usize)(x);

        let node_disps = self.get_nodes_displacement();
        let u = n0 * node_disps[0] + n1 * node_disps[1];
        [u]
    }

    /// Return a 4x4 element stiffness matrix, elements are Dtype
    fn calc_k(&self) -> [[Dtype; 4]; 4] {
        println!(
            "\n>>> Calculating Rod2D2N(#{})'s stiffness matrix k{} ......",
            self.id, self.id
        );
        let l = self.length();
        let ee = self.material.0;
        let tmat = self.trans_mat();
        let area = self.cross_sectional_area;
        let base = [[1.0, -1.0], [-1.0, 1.0]];
        type SM2 = SMatrix<Dtype, 2, 2>;
        let local_k_mat = (SM2::from(base)) * (ee * area / l);
        (tmat.transpose() * local_k_mat * tmat).into()
    }

    /// Get element's strain vector, in 1d it's a scale
    fn calc_strain(&self) -> [Dtype; 1] {
        let l = self.length();
        let tmat = self.trans_mat();
        let disp = self.get_nodes_displacement();
        let dvec = SMatrix::<Dtype, 4, 1>::from(disp);
        let bmat = SMatrix::<Dtype, 1, 2>::from([-1.0 / l, 1.0 / l]);
        let strain: [Dtype; 1] = (bmat * tmat * dvec).into();
        strain
    }

    /// Get element's stress vector, in 1d it's a scale
    fn calc_stress(&self) -> [Dtype; 1] {
        let ee = self.material.0;
        let strain = self.calc_strain()[0];
        [ee * strain]
    }

    /// Output element's calculation result
    pub fn calc_result_info(&self, n_exp: Dtype) -> String {
        format!("\n-----------------------------------------------------------------------------\nElem_Rod2d2N:\n\tId:\t{}\n\tSection Area: {:-12.6}\n\tMats: {:-12.6} (Young's modulus)\n\t      {:-12.6} (Poisson's ratio)\n\tNodes:{}{}\n\tStrain:\n\t\t{:-12.6?}\n\tStress:\n\t\t{:-12.6?}\n\n\tStiffness Matrix K{} =  (*10^{})\n{}",
            self.id,
            self.cross_sectional_area,
            self.material.0,
            self.material.1,
            self.nodes[0],
            self.nodes[1],
            self.calc_strain(),
            self.calc_stress(),
            self.id(),
            n_exp,
            self.k_string(n_exp),
        )
    }

    /// Print element's strain value
    pub fn print_strain(&self) {
        let strain = self.calc_strain();
        println!(
            "\nelem[{}] strain:\n\tE_xx = {:-12.6}\n\tE_yy = {:-9.6}\n\tE_xy = {:-9.6}",
            self.id, strain[0], 0.0, 0.0
        );
    }

    /// Print element's stress value
    pub fn print_stress(&self) {
        let stress = self.calc_stress();
        println!("\nelem[{}] stress:\n\tS_xx = {:-12.6}", self.id, stress[0]);
    }
}

/// Implement zhm::K trait for Rod1D2N element
impl<'rod2d2n> K for Rod2D2N<'rod2d2n> {
    /// Cache stiffness matrix for rod element
    fn k(&mut self) -> &CompressedMatrix {
        if self.k_matrix.is_none() {
            self.k_matrix.get_or_insert(compress_matrix(self.calc_k()))
        } else {
            self.k_matrix.as_ref().unwrap()
        }
    }

    /// Print Rod1D2N element's stiffness matrix
    fn k_printer(&self, n_exp: Dtype) {
        if self.k_matrix.is_none() {
            panic!(
                "!!! Rod1D2N#{}'s k mat is empty! call k() to calc it.",
                self.id
            );
        }
        print_2darr("\nRod2d2N k", self.id(), &self.calc_k(), n_exp);
    }

    /// Return Rod1D2N elem's stiffness matrix's format string
    fn k_string(&self, n_exp: Dtype) -> String {
        let mut k_matrix = String::new();
        let elem_stiffness_mat: [[Dtype; 4]; 4] = self.calc_k();
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
                    elem_stiffness_mat[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
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

    /// Get the strain at (x,y) inside the element, in linear rod elem, strain is a const
    fn strain_at_intpt(&mut self) -> Data {
        if self.strain.is_none() {
            self.strain.get_or_insert(self.calc_strain());
        }
        Data::Dim1(vec![self.strain.unwrap()])
    }

    /// Get the stress at (x,y) insDatae the element, in linear rod elem, stress is a const
    fn stress_at_intpt(&mut self) -> Data {
        if self.stress.is_none() {
            self.stress.get_or_insert(self.calc_stress());
        }
        Data::Dim1(vec![self.stress.unwrap()])
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

impl fmt::Display for Rod2D2N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nElem_Rod2d2N:\n\tId:\t{}\n\tArea: {:-12.4}\n\tMats: {:-12.4} (Young's modulus)\n\t      {:-12.4} (Poisson's ratio)\n\tNodes:{}{}",
            self.id,
            self.cross_sectional_area,
            self.material.0,
            self.material.1,
            self.nodes[0],
            self.nodes[1],
        )
    }
}
