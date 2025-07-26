use crate::dtty::{
    basic::{Dtype, Jacobian2D},
    matrix::CompressedMatrix,
};
use crate::node::Node2D;
use crate::port::K;
use crate::tool::{compress_matrix, print_2darr};
use na::SMatrix;
use std::fmt::{self, Write};

/// Three-node triangular element in two-dimensional plane
/// The stress and strain in the CST element are constants
pub struct Tri2D3N<'tri2d3n> {
    pub id: usize,
    pub thick: Dtype,
    pub nodes: [&'tri2d3n Node2D; 3],
    pub strain: Option<[Dtype; 3]>,
    pub stress: Option<[Dtype; 3]>,
    pub k_matrix: Option<CompressedMatrix>,
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
            strain: None,
            stress: None,
            k_matrix: None,
            material,
        }
    }

    /// Set Tri2D3N element's material parameter
    pub fn set_material(&mut self, material_args: &'tri2d3n (Dtype, Dtype)) {
        self.material = material_args;
    }

    /// Get Tri2D3N element's area value
    pub fn area(&self) -> Dtype {
        let x = self.get_nodes_xcoords();
        let y = self.get_nodes_ycoords();
        let dx_21 = x[1] - x[0];
        let dx_31 = x[2] - x[0];
        let dy_21 = y[1] - y[0];
        let dy_31 = y[2] - y[0];
        0.5 * (dx_21 * dy_31 - dx_31 * dy_21).abs()
    }

    /// Get the x-coords of nodes in Tri2D3N element
    pub fn get_nodes_xcoords(&self) -> [Dtype; 3] {
        let mut x_list = [0.0; 3];
        for (idx, &node) in self.nodes.iter().enumerate() {
            x_list[idx] = node.coords[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in Tri2D3N element
    pub fn get_nodes_ycoords(&self) -> [Dtype; 3] {
        let mut y_list = [0.0; 3];
        for (idx, &node) in self.nodes.iter().enumerate() {
            y_list[idx] = node.coords[1];
        }
        y_list
    }

    /// Transform area-coords [L1, L2, L3] into (x, y) coords
    pub fn coords_transform_ltox(&self, area_coords: [Dtype; 2]) -> [Dtype; 2] {
        let xs = self.get_nodes_xcoords();
        let ys = self.get_nodes_ycoords();
        let l3 = 1.0 - area_coords[0] - area_coords[1];
        let x = xs[0] * area_coords[0] + xs[1] * area_coords[1] + xs[2] * l3;
        let y = ys[0] * area_coords[0] + ys[1] * area_coords[1] + ys[2] * l3;
        [x, y]
    }

    /// Get nodes' disps vector in tri element
    pub fn get_nodes_displacement(&self) -> [Dtype; 6] {
        let mut displacement = [0.0; 6];
        for (idx, &node) in self.nodes.iter().enumerate() {
            displacement[2 * idx] = node.displs.borrow()[0];
            displacement[2 * idx + 1] = node.displs.borrow()[1];
        }
        displacement
    }

    /// Get any point's disps vector in tri element
    pub fn get_point_displacement(&self, area_coords: [Dtype; 2]) -> [Dtype; 2] {
        let xy = self.coords_transform_ltox(area_coords);
        let x = xy[0];
        let y = xy[1];
        let n0 = self.shape_func_xy(0usize)(x, y);
        let n1 = self.shape_func_xy(1usize)(x, y);
        let n2 = self.shape_func_xy(2usize)(x, y);

        let disps = self.get_nodes_displacement();
        let u = n0 * disps[0] + n1 * disps[2] + n2 * disps[4];
        let v = n0 * disps[1] + n1 * disps[3] + n2 * disps[5];
        [u, v]
    }

    /// Get nodes's force vector in tri element
    pub fn get_nodes_force(&self) -> [Dtype; 6] {
        let mut force = [0.0; 6];
        for (idx, &node) in self.nodes.iter().enumerate() {
            force[2 * idx] = node.forces.borrow()[0];
            force[2 * idx + 1] = node.forces.borrow()[1];
        }
        force
    }

    /// Coefficient of shape function N(x, y)
    /// u(x, y) = A0 + A1 * x + A2 * y
    /// v(x, y) = B0 + B1 * x + B2 * y
    ///
    /// a,b,c通过计算下式首行的代数余子式获得：
    /// | 1  X1  Y1 |
    /// | 1  X2  Y2 |
    /// | 1  X3  Y3 |
    /// 具体的表达式为：
    /// Ai =  (Xj * Ym - Xm * Yj)
    /// Bi =  (Yj      - Ym     )  (i -> j -> m)
    /// Ci = -(Xj      - Xm     )
    ///
    /// Ni = (1/2A) * (ai + bi * x + ci * y) (i -> j -> m)
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
    /// Bi =  Yj - Ym  (i -> j -> m)
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
    /// Ci = -(Xj - Xm)
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
    /// The shape mat of Tri2D3N:
    /// |N1  0   N2  0   N3  0 |
    /// |0   N1  0   N2  0   N3|
    /// Ni = (1/2A) * (ai + bi * x + ci * y) (i -> j -> m)
    /// 这里的A是三角形单元的面积
    fn shape_func_xy(&self, ith: usize) -> impl Fn(Dtype, Dtype) -> Dtype {
        let area = self.area();
        let a = self.coef_a(ith);
        let b = self.coef_b(ith);
        let c = self.coef_c(ith);
        move |x: Dtype, y: Dtype| 0.5 * (a + x * b + y * c) / area
    }

    /// Derivatives of shape functions with respect to physical coordinates
    /// 形函数对物理坐标的导数
    /// dN/dx = (dN/ds * ds/dx) + (dN/dt * dt/dx)
    /// dN/dy = (dN/ds * ds/dy) + (dN/dt * dt/dy)
    /// | dN/ds |        | dx/ds  dy/ds | | dN/dx |          | dN/dx |
    /// |       |    =   |              |*|       |   =  J * |       |
    /// | dN/dt |        | dx/dt  dy/dt | | dN/dy |2x3       | dN/dy |2x3
    ///
    /// | dN/dx |                 | dN/ds |
    /// |       |    =   J^(-1) * |       |
    /// | dN/dy |2x3              | dN/dt |2x3
    #[inline]
    fn diff_shape_mat_xy(&self) -> SMatrix<Dtype, 2, 3> {
        let diff_shape_mat_st =
            SMatrix::<Dtype, 2, 3>::from([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]]);
        let jecobian = self.jacobian();
        jecobian.try_inverse().unwrap() * diff_shape_mat_st
    }

    /// B(s,t) for Tri2D3N which is a 3x6 mat
    #[inline]
    fn geometry_mat_xy(&self) -> SMatrix<Dtype, 3, 6> {
        let dn_xy = self.diff_shape_mat_xy();
        SMatrix::<Dtype, 3, 6>::from([
            [dn_xy[(0, 0)], 0.0, dn_xy[(1, 0)]],
            [0.0, dn_xy[(1, 0)], dn_xy[(0, 0)]],
            [dn_xy[(0, 1)], 0.0, dn_xy[(1, 1)]],
            [0.0, dn_xy[(1, 1)], dn_xy[(0, 1)]],
            [dn_xy[(0, 2)], 0.0, dn_xy[(1, 2)]],
            [0.0, dn_xy[(1, 2)], dn_xy[(0, 2)]],
        ])
    }

    /// Calculate the Jacobian matrix of tri2d3n element
    /// J = [[dx/ds, dy/ds]
    ///      [dx/dt, dy/dt]]
    #[inline]
    fn jacobian(&self) -> Jacobian2D {
        let xs: [Dtype; 3] = self.get_nodes_xcoords();
        let ys: [Dtype; 3] = self.get_nodes_ycoords();
        let dx_21 = xs[1] - xs[0];
        let dx_31 = xs[2] - xs[0];
        let dy_21 = ys[1] - ys[0];
        let dy_31 = ys[2] - ys[0];
        Jacobian2D::from([[dx_21, dy_21], [dx_31, dy_31]]).transpose()
    }

    /// Return a 6x6 matrix, elements are Dtype
    /// K(x, y) = thick * INT[ B(x, y)' * D * B(x, y) ] dA
    ///         = thick * INT[ B(s, t)' * D * B(s, t) * det(J) ] dsdt
    ///         = thick * B(s, t)' * D * B(s, t) * det(J) * Area
    ///         = thick * B(s, t)' * D * B(s, t) * det(J) * 0.5 //标准三角形面积为0.5
    fn calc_k(&self) -> [[Dtype; 6]; 6] {
        //println!(
        //    "\n>>> Calculating Tri2D3N(#{})'s local stiffness matrix k{} ......",
        //    self.id, self.id
        //);
        let (ee, nu) = *self.material;
        let elasticity_mat = (ee / (1.0 - nu * nu))
            * (SMatrix::<Dtype, 3, 3>::from([
                [1.0, nu, 0.0],
                [nu, 1.0, 0.0],
                [0.0, 0.0, 0.5 * (1.0 - nu)],
            ])
            .transpose());
        let det_j = self.jacobian().determinant();

        let b_mat = self.geometry_mat_xy();
        let core = b_mat.transpose() * elasticity_mat * b_mat * det_j;
        let stiffness_matrix: [[Dtype; 6]; 6] = (0.5 * self.thick * core).into();
        stiffness_matrix
    }

    /// Get element's strain vector, the strain in CST elem is a const
    pub fn calc_strain(&self) -> [Dtype; 3] {
        let b_mat = self.geometry_mat_xy();
        let elem_nodes_disps = SMatrix::<Dtype, 6, 1>::from(self.get_nodes_displacement());
        (b_mat * elem_nodes_disps).into()
    }

    /// Get element's stress vector, the stress in CST elem is a const
    pub fn calc_stress(&self) -> [Dtype; 3] {
        let (ee, nu) = *self.material;
        let elasticity_mat = SMatrix::<Dtype, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        let epsilon = SMatrix::<Dtype, 3, 1>::from(self.calc_strain());
        (elasticity_mat * epsilon).into()
    }

    /// Output element's calculation result
    pub fn calc_result_info(&self, n_exp: Dtype) -> String {
        format!("\n-----------------------------------------------------------------------------\nElem_Tri2D3N:\n\tId:\t{}\n\tArea: {:-12.6}\n\tMats: {:-12.6} (Young's modulus)\n\t      {:-12.6} (Poisson's ratio)\n\tNodes:{}{}{}\n\tStrain:\n\t\t{:-12.6?}\n\tStress:\n\t\t{:-12.6?}\n\n\tStiffness Matrix K{} =  (*10^{})\n{}",
            self.id,
            self.area(),
            self.material.0,
            self.material.1,
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
            self.calc_strain(),
            self.calc_stress(),
            self.id(),
            n_exp,
            self.k_string(n_exp),
        )
    }

    /// Print element's strain value
    pub fn print_strain(&self) {
        let strain = self.strain.unwrap();
        println!(
            "\nelem[{}] strain:\n\tE_xx = {:-16.6}\n\tE_yy = {:-16.6}\n\tE_xy = {:-16.6}",
            self.id, strain[0], strain[1], strain[2]
        );
    }

    /// Print element's stress value
    pub fn print_stress(&self) {
        let stress = self.stress.unwrap();
        println!(
            "\nelem[{}] stress:\n\tS_xx = {:-16.6}\n\tS_yy = {:-16.6}\n\tS_xy = {:-16.6}",
            self.id, stress[0], stress[1], stress[2]
        );
    }
}

/// Implement zhm::K trait for triangle element
impl<'tri2d3n> K for Tri2D3N<'tri2d3n> {
    /// Cache stiffness matrix for triangle element
    fn k(&mut self) -> &CompressedMatrix {
        if self.k_matrix.is_none() {
            self.k_matrix.get_or_insert(compress_matrix(self.calc_k()))
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
        print_2darr("\nTri2D3N k", self.id(), &self.calc_k(), n_exp);
    }

    /// Return triangle elem's stiffness matrix's format string
    fn k_string(&self, n_exp: Dtype) -> String {
        let mut k_matrix = String::new();
        let elem_stiness_mat: [[Dtype; 6]; 6] = self.calc_k();
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
                    elem_stiness_mat[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
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

    /// Get the strain at integration point
    fn strain_at_intpt(&mut self) -> Vec<Vec<Dtype>> {
        if self.strain.is_none() {
            self.strain = Some(self.calc_strain());
        }
        vec![self.strain.unwrap().to_vec()]
    }

    /// Get the stress at integratiData point
    fn stress_at_intpt(&mut self) -> Vec<Vec<Dtype>> {
        if self.stress.is_none() {
            let cst_stress = self.calc_stress();
            self.stress = Some(cst_stress);
        }
        vec![self.stress.unwrap().to_vec()]
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

impl fmt::Display for Tri2D3N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let strain = self.strain.unwrap();
        let stress = self.stress.unwrap();
        write!(
            f,
"\n-----------------------------------------------------------------------------\nElem_Tri2D3N:\n\tId:\t{}\n\tArea: {:-12.6}\n\tMats: {:-12.6} (Young's modulus)\n\t      {:-12.6} (Poisson's ratio)\n\tNodes:{}{}{}\n\tStrain:\n\t\t{:-12.6?}\n\tStress:\n\t\t{:-12.6?}\n\tStiffness Matrix K{} =  (*10^0)\n{}",
            self.id,
            self.area(),
            self.material.0,
            self.material.1,
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
            strain,
            stress,
            self.id(),
            self.k_string(0.0),
        )
    }
}

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
    pub fn get_point_displacement(&self, point_coord_local: [Dtype; 3]) -> [Dtype; 2] {
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
