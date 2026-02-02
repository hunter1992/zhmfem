use crate::dtty::{basic::Dtype, matrix::CompressedMatrixSKS};
use crate::node::Node2D;
use crate::port::SData;
use crate::port::StaticStiffness;
use crate::tool::{compress_symmetry_matrix_sks, print_2darr};
use na::SMatrix;
use std::fmt::{self, Write};

pub struct Rod2D2N<'rod2d2n> {
    pub id: usize,
    pub cross_section_area: Dtype,
    pub nodes: [&'rod2d2n Node2D; 2],
    pub k_matrix: Option<CompressedMatrixSKS>,
    pub material: [Dtype; 2],
}

impl<'rod2d2n> Rod2D2N<'rod2d2n> {
    /// Generate a 2D Rod2D2N element
    pub fn new(
        id: usize,
        cross_section_area: Dtype,
        nodes: [&'rod2d2n Node2D; 2],
        material: [Dtype; 2],
    ) -> Self {
        Rod2D2N {
            id,
            cross_section_area,
            nodes,
            k_matrix: None,
            material,
        }
    }

    /// Get element's id number
    fn id(&self) -> usize {
        self.id
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

    /// Get Rod2D2N element's length value
    fn length(&self) -> Dtype {
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

    /// Set element material_args
    pub fn set_material(&mut self, material_args: [Dtype; 2]) {
        self.material = material_args;
    }

    /// Get the x-coords of two nodes in Rod1D2N element
    pub fn get_nodes_xcoords(&self) -> [Dtype; 2] {
        let mut x_list: [Dtype; 2] = [0.0; 2];
        for idx in 0..2 {
            x_list[idx] = self.nodes[idx].coords[0];
        }
        x_list
    }

    /// Get the y-coords of two nodes in Rod1D2N element
    pub fn get_nodes_ycoords(&self) -> [Dtype; 2] {
        let mut y_list: [Dtype; 2] = [0.0; 2];
        for idx in 0..2 {
            y_list[idx] = self.nodes[idx].coords[1];
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
        self.calc_stress()[0] * self.cross_section_area
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
        let ee = self.material[0];
        let tmat = self.trans_mat();
        let area = self.cross_section_area;
        let base = [[1.0, -1.0], [-1.0, 1.0]];
        type SM2 = SMatrix<Dtype, 2, 2>;
        let local_k_mat = (SM2::from(base)) * (ee * area / l);
        (tmat.transpose() * local_k_mat * tmat).into()
    }

    /// Get element's strain vector
    fn calc_strain(&self) -> [Dtype; 1] {
        let l = self.length();
        let tmat = self.trans_mat();
        let disp = self.get_nodes_displacement();
        let dvec = SMatrix::<Dtype, 4, 1>::from(disp);
        let bmat = SMatrix::<Dtype, 1, 2>::from([-1.0 / l, 1.0 / l]);
        let strain: [Dtype; 1] = (bmat * tmat * dvec).into();
        strain
    }

    /// Get element's stress vector
    fn calc_stress(&self) -> [Dtype; 1] {
        let ee = self.material[0];
        let strain = self.calc_strain()[0];
        [ee * strain]
    }

    /// Output element's calculation result
    pub fn calc_result_info(&self, n_exp: Dtype) -> String {
        format!("\n-----------------------------------------------------------------------------\nElem_Rod2d2N:\n\tId:\t{}\n\tSection Area: {:-12.6}\n\tMats: {:-12.6} (Young's modulus)\n\t      {:-12.6} (Poisson's ratio)\n\tNodes:{}{}\n\tStrain:\n\t\t{:-12.6?}\n\tStress:\n\t\t{:-12.6?}\n\n\tStiffness Matrix K{} =  (*10^{})\n{}",
            self.id,
            self.cross_section_area,
            self.material[0],
            self.material[1],
            self.nodes[0],
            self.nodes[1],
            self.calc_strain(),
            self.calc_stress(),
            self.id(),
            n_exp,
            self.k_string(n_exp),
        )
    }
}

impl<'rod2d2n> StaticStiffness for Rod2D2N<'rod2d2n> {
    /// Cache stiffness matrix for rod element
    fn k(&mut self) -> &CompressedMatrixSKS {
        if self.k_matrix.is_none() {
            self.k_matrix
                .get_or_insert(compress_symmetry_matrix_sks(&self.calc_k()))
        } else {
            self.k_matrix.as_ref().unwrap()
        }
    }

    /// Print Rod2D2N element's stiffness matrix
    fn k_printr(&self, n_exp: Dtype) {
        let k_mat: [[Dtype; 4]; 4] = self.calc_k();
        print_2darr("\nRod2D2N k", self.id, &k_mat, n_exp)
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
}

impl<'rod2d2n> SData for Rod2D2N<'rod2d2n> {
    /// Get rod1d2n element's id
    fn nodes_ids(&self) -> Vec<usize> {
        let mut id: Vec<usize> = Vec::with_capacity(2);
        for idx in 0..2 {
            id.push(self.nodes[idx].id);
        }
        id
    }

    /// Get Rod1D2N element length
    fn elem_size(&self) -> Dtype {
        self.length()
    }

    /// Get the strain at x inside the element, in linear rod elem, strain is a const
    fn strain_at_nodes(&mut self) -> Vec<Dtype> {
        self.calc_strain().to_vec()
    }

    /// Get the stress at x inside the element, in linear rod elem, stress is a const
    fn stress_at_nodes(&mut self) -> Vec<Dtype> {
        self.calc_stress().to_vec()
    }
}

impl fmt::Display for Rod2D2N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nElem_Rod2d2N:\n\tId:\t{}\n\tArea: {:-12.4}\n\tMats: {:-12.4} (Young's modulus)\n\t      {:-12.4} (Poisson's ratio)\n\tNodes:{}{}",
            self.id,
            self.cross_section_area,
            self.material[0],
            self.material[1],
            self.nodes[0],
            self.nodes[1],
        )
    }
}
