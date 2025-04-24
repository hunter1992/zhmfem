use crate::{compress_matrix, node::Node1D, CompressedMatrix, Data, Dtype, K};
use na::SMatrix;
use std::fmt::{self, Write};

pub struct Rod1D2N<'rod1d2n> {
    pub id: usize,
    pub cross_sectional_area: Dtype,
    pub nodes: [&'rod1d2n Node1D; 2],
    pub strain: Option<[Dtype; 1]>,
    pub stress: Option<[Dtype; 1]>,
    pub k_matrix: Option<CompressedMatrix>,
    pub material: &'rod1d2n (Dtype, Dtype),
}

impl<'rod1d2n> Rod1D2N<'rod1d2n> {
    /// Generate a 1D Rod1D2N element
    pub fn new(
        id: usize,
        cross_sectional_area: Dtype,
        nodes: [&'rod1d2n Node1D; 2],
        material: &'rod1d2n (Dtype, Dtype),
    ) -> Self {
        Rod1D2N {
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
    pub fn set_material(&mut self, material_args: &'rod1d2n (Dtype, Dtype)) {
        self.material = material_args;
    }

    /// Get Rod1D2N element length
    pub fn length(&self) -> Dtype {
        let x = self.get_nodes_xcoords();
        (x[0] - x[1]).abs()
    }

    /// Get the x-coords of nodes in Rod1D2N element
    pub fn get_nodes_xcoords(&self) -> [Dtype; 2] {
        let mut x_list = [0.0; 2];
        for i in 0..2 {
            x_list[i] = self.nodes[i].coords[0];
        }
        x_list
    }

    /// Get nodes' disps vector in Rod1D2N element
    pub fn get_nodes_displacement(&self) -> [Dtype; 2] {
        let mut disps = [0.0; 2];
        for idx in 0..2 {
            disps[idx] = self.nodes[idx].displs.borrow()[0];
        }
        disps
    }

    /// Get nodes's force vector in Rod1D2N element
    pub fn get_nodes_force(&self) -> [Dtype; 2] {
        let mut forces = [0.0; 2];
        for idx in 0..2 {
            forces[idx] = self.nodes[idx].forces.borrow()[0];
        }
        forces
    }

    /// Get axial force in Rod1D2N element
    pub fn axial_force(&self) -> Dtype {
        self.calc_stress()[0] * self.cross_sectional_area
    }

    /// Get shape matrix element N_i
    fn shape_mat_i(&self, ith: usize) -> impl Fn(Dtype) -> Dtype {
        /* The shape mat of rod elem:
         * [N1 N2]  which is a 1x2 mat
         * N1 = 1 - x/L = 1 - epsilon
         * N2 = x/L     = epsilon      */
        let a: [Dtype; 2] = [1.0, 0.0];
        let b: [Dtype; 2] = [-1.0, 1.0];
        let length = self.length();
        // a[0] + b[0]*epsilon构造出N1,a[1] + b[1]*epsilon构造出N2
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

    /// Return a 2x2 element stiffness matrix, elements are Dtype
    fn calc_k(&self) -> [[Dtype; 2]; 2] {
        println!(
            "\n>>> Calculating Rod1D2N(#{})'s stiffness matrix k{} ......",
            self.id, self.id
        );
        let ee = self.material.0;
        let stiffness_matrix: [[Dtype; 2]; 2] =
            (SMatrix::<Dtype, 2, 2>::from([[1.0, -1.0], [-1.0, 1.0]])
                * (ee * self.cross_sectional_area / self.length()))
            .into();
        stiffness_matrix
    }

    /// Get element's strain vector, in 1d it's a scale
    fn calc_strain(&self) -> [Dtype; 1] {
        let unit: Dtype = 1.0 / self.length();
        let node_disps = self.get_nodes_displacement();
        let strain: Dtype = -unit * node_disps[0] + unit * node_disps[1];
        [strain]
    }

    /// Get element's stress vector, in 1d it's a scale
    fn calc_stress(&self) -> [Dtype; 1] {
        let ee = self.material.0;
        let stress: Dtype = ee * self.calc_strain()[0];
        [stress]
    }

    /// Print element's strain value
    pub fn print_strain(&self) {
        let strain = self.calc_strain();
        println!(
            "\nelem[{}] strain:\n\tE_xx = {:-12.6}\n\t",
            self.id, strain[0]
        );
    }

    /// Print element's stress value
    pub fn print_stress(&self) {
        let stress = self.calc_stress();
        println!(
            "\nelem[{}] stress:\n\tS_xx = {:-12.6}\n\t",
            self.id, stress[0]
        );
    }
}

/// Implement zhm::K trait for Rod1D2N element
impl<'rod1d2n> K for Rod1D2N<'rod1d2n> {
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

        print!("\nRod1D2N k{} =  (* 10^{})\n[", self.id, n_exp as u8);
        for row in 0..2 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..2 {
                print!(
                    " {:>-12.6} ",
                    self.calc_k()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                );
            }
            if row == 1 {
                println!("]]");
            } else {
                println!("]");
            }
        }
        print!("\n");
    }

    /// Return Rod1D2N elem's stiffness matrix's format string
    fn k_string(&self, n_exp: Dtype) -> String {
        let mut k_matrix = String::new();
        let elem_stiffness_mat = self.calc_k();
        for row in 0..2 {
            if row == 0 {
                write!(k_matrix, "[[").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, " [").expect("!!! Write tri k_mat failed!");
            }
            for col in 0..2 {
                write!(
                    k_matrix,
                    " {:>-12.6} ",
                    elem_stiffness_mat[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                )
                .expect("!!! Write tri k_mat failed!");
            }
            if row == 1 {
                write!(k_matrix, "]]").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, "]\n").expect("!!! Write tri k_mat failed!");
            }
        }
        k_matrix
    }

    /// Get the strain at (x,y) inside the element, in linear rod elem, strain is a const
    fn strain_intpt(&mut self) -> Data {
        if self.strain.is_none() {
            self.strain.get_or_insert(self.calc_strain());
        }
        Data::Dim1(vec![self.strain.unwrap()])
    }

    /// Get the stress at (x,y) insDatae the element, in linear rod elem, stress is a const
    fn stress_intpt(&mut self) -> Data {
        if self.stress.is_none() {
            self.stress.get_or_insert(self.calc_stress());
        }
        Data::Dim1(vec![self.stress.unwrap()])
    }

    /// Get element's info string
    fn info(&self, n_exp: Dtype) -> String {
        format!("\n-----------------------------------------------------------------------------\nElem_Rod1d2N:\n\tId:\t{}\n\tArea: {:-12.4} (cross sectional area)\n\tLen : {:-12.4}\n\tMats: {:-12.4} (Young's modulus)\n\t      {:-12.4} (Poisson's ratio)\n\tNodes:{}{}\n\tStrain:\n\t\t{:-12.6?}\n\n\tStress:\n\t\t{:-12.6?}\n\n\tStiffness Matrix K{} =  (* 10^{})\n{}",
            self.id,
            self.cross_sectional_area,
            self.length(),
            self.material.0,
            self.material.1,
            self.nodes[0],
            self.nodes[1],
            self.calc_strain(),
            self.calc_stress(),
            self.id(),
            n_exp,
            self.k_string(n_exp)
        )
    }

    /// Get element's id number
    fn id(&self) -> usize {
        self.id
    }
}

impl fmt::Display for Rod1D2N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nElem_Rod1d2N:\n\tId:\t{}\n\tArea: {:-12.4}\n\tMats: {:-12.4} (Young's modulus)\n\t      {:-12.4} (Poisson's ratio)\n\tNodes:{}{}",
            self.id,
            self.cross_sectional_area,
            self.material.0,
            self.material.1,
            self.nodes[0],
            self.nodes[1],
        )
    }
}
