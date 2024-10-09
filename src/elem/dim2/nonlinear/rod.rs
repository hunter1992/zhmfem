use crate::{node::Node2D, Dtype, K};
use na::*;
use std::fmt::Write;

pub struct Rod2D2NNL<'rod> {
    pub id: usize,
    pub sec_area: Dtype,
    pub nodes: [&'rod Node2D; 2],
    pub tangent_kmatrix: Option<[[Dtype; 4]; 4]>,
}

#[allow(non_snake_case)]
impl<'rod> Rod2D2NNL<'rod> {
    /// Generate a 2D Rod2D2N element
    pub fn new(id: usize, sec_area: Dtype, nodes: [&Node2D; 2]) -> Rod2D2NNL {
        Rod2D2NNL {
            id,
            sec_area,
            nodes,
            tangent_kmatrix: None,
        }
    }

    /// Get Node1's X coordinate under initial configuration
    pub fn X1(&self) -> Dtype {
        let X1: Dtype = self.nodes[0].coords[0];
        X1
    }

    /// Get Node1's x-coord displacement
    pub fn disp1_x(&self) -> Dtype {
        let disp1_x = self.nodes[0].displs[0];
        disp1_x
    }

    /// Get Node1's x coordinate under current configuration
    pub fn x1(&self) -> Dtype {
        let x1: Dtype = self.X1() + self.disp1_x();
        x1
    }

    /// Get Node1's Y coordinate under initial configuration
    pub fn Y1(&self) -> Dtype {
        let Y1: Dtype = self.nodes[0].coords[1];
        Y1
    }

    /// Get Node1's y-coord displacement
    pub fn disp1_y(&self) -> Dtype {
        let disp1_y = self.nodes[0].displs[1];
        disp1_y
    }

    /// Get Node1's y coordinate under current configuration
    pub fn y1(&self) -> Dtype {
        let y1: Dtype = self.Y1() + self.disp1_y();
        y1
    }

    /// Get Node2's X coordinate under initial configuration
    pub fn X2(&self) -> Dtype {
        let X2: Dtype = self.nodes[1].coords[0];
        X2
    }

    /// Get Node2's x-coord displacement
    pub fn disp2_x(&self) -> Dtype {
        let disp2_x = self.nodes[1].displs[0];
        disp2_x
    }

    /// Get Node2's x coordinate under current configuration
    pub fn x2(&self) -> Dtype {
        let x2: Dtype = self.X2() + self.disp2_x();
        x2
    }

    /// Get Node2's Y coordinate under initial configuration
    pub fn Y2(&self) -> Dtype {
        let Y2: Dtype = self.nodes[1].coords[1];
        Y2
    }

    /// Get Node2's y-coord displacement
    pub fn disp2_y(&self) -> Dtype {
        let disp2_y = self.nodes[1].displs[1];
        disp2_y
    }

    /// Get Node2's y coordinate under current configuration
    pub fn y2(&self) -> Dtype {
        let y2: Dtype = self.Y2() + self.disp2_y();
        y2
    }

    /// Get delta X under initial configuration
    pub fn dX(&self) -> Dtype {
        let dX: Dtype = self.X2() - self.X1();
        dX
    }

    /// Get delta x under current configuration
    pub fn dx(&self) -> Dtype {
        let dx: Dtype = self.x2() - self.x1();
        dx
    }

    /// Get delta Y under initial configuration
    pub fn dY(&self) -> Dtype {
        let dY: Dtype = self.Y2() - self.Y1();
        dY
    }

    /// Get delta y under current configuration
    pub fn dy(&self) -> Dtype {
        let dy: Dtype = self.y2() - self.y1();
        dy
    }

    /// Get rod element length under initial configuration
    pub fn length_init(&self) -> Dtype {
        let dX: Dtype = self.dX();
        let dY: Dtype = self.dY();
        let l_init: Dtype = (dX * dX + dY * dY).sqrt();
        l_init
    }

    /// Get rod element length under current configuration
    pub fn length_current(&self) -> Dtype {
        let dx: Dtype = self.dx();
        let dy: Dtype = self.dy();
        let l_current: Dtype = (dx * dx + dy * dy).sqrt();
        l_current
    }

    /// Get sin(Alpha) under origin configuration,
    /// Alpha is the angle between X coord and rod length direction
    pub fn sin_Alpha(&self) -> Dtype {
        self.dY() / self.length_init()
    }

    /// Get cos(Alpha) under origin configuration,
    /// Alpha is the angle between X coord and rod length direction
    pub fn cos_Alpha(&self) -> Dtype {
        self.dX() / self.length_init()
    }

    /// Get sin(alpha) under origin configuration,
    /// alpha is the angle between x coord and rod length direction
    pub fn sin_alpha(&self) -> Dtype {
        self.dy() / self.length_current()
    }

    /// Get cos(alpha) under origin configuration,
    /// alpha is the angle between x coord and rod length direction
    pub fn cos_alpha(&self) -> Dtype {
        self.dx() / self.length_current()
    }

    /// Variation of strain to displacement is B mat,
    /// B matrix  = d(strain)/d(disp)
    ///           = (1/L0) * [-ax, -ay, ax, ay], which is 1x4 mat
    /// and   ax  =  (x2 - x1) / L0,
    ///       ay  =  (y2 - y1) / L0
    pub fn strain_variation_mat(&self) -> [Dtype; 4] {
        let ax = (self.x2() - self.x1()) / self.length_init();
        let ay = (self.y2() - self.y1()) / self.length_init();
        [-ax, -ay, ax, ay]
    }

    /// Get Green-Lagrange strain of rod element
    pub fn green_langrage_strain(&self) -> [Dtype; 3] {
        let l_init: Dtype = self.length_init();
        let l_current: Dtype = self.length_current();
        let strain_gl: Dtype =
            (l_current * l_current - l_init * l_init) / (2.0 as Dtype * l_init * l_init);
        [strain_gl, 0.0, 0.0]
    }

    /// Get the Second Piola-Kirchhoff Stress
    pub fn PKII_stress(&self, material_args: (Dtype, Dtype)) -> [Dtype; 3] {
        let pk2_stess = material_args.0 * self.green_langrage_strain()[0];
        [pk2_stess, 0.0, 0.0]
    }

    /// Calculate Rod element's tangent stiffness matrix Kt under nonlinear analysis.
    /// Return a 4x4 matrix, elements are Dtype
    fn tangent_stiffness_mat(&self, material_args: (Dtype, Dtype)) -> [[Dtype; 4]; 4] {
        println!(
            "\n>>> Calculating Rod2D2NNL(#{})'s tangent stiffness matrix Kt{}",
            self.id, self.id
        );

        let km_mat =
            SMatrix::<Dtype, 4, 4>::from(self.material_stiffness_mat(material_args)).transpose();
        let kg_mat =
            SMatrix::<Dtype, 4, 4>::from(self.geometric_stiffness_mat(material_args)).transpose();
        let tangent_kmat: [[Dtype; 4]; 4] = (km_mat + kg_mat).into();
        tangent_kmat
    }

    /// Get material stiffness matrix k_m of rod element under nonlinear analysis
    pub fn material_stiffness_mat(&self, material_args: (Dtype, Dtype)) -> [[Dtype; 4]; 4] {
        let ee = material_args.0;
        let L0 = self.length_init();
        let strain_mat = SMatrix::<Dtype, 1, 4>::from(self.strain_variation_mat());
        let material_kmat: [[Dtype; 4]; 4] =
            (ee / L0 * self.sec_area * (strain_mat.transpose() * strain_mat)).into();
        material_kmat
    }

    /// Get geometry stiffness matrix k_g of rod element under nonlinear analysis
    pub fn geometric_stiffness_mat(&self, material_args: (Dtype, Dtype)) -> [[Dtype; 4]; 4] {
        let ee = material_args.0;
        let b_mat = self.green_langrage_strain()[0];
        let arg: Dtype = ee * self.sec_area * b_mat / self.length_init();
        let base_mat = SMatrix::<Dtype, 4, 4>::from([
            [1.0, 0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0, -1.0],
            [-1.0, 0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0, 1.0],
        ]);
        let geometric_kmat: [[Dtype; 4]; 4] = (arg * base_mat).into();
        geometric_kmat
    }

    /*
    /// node's inner force writer
    /// inner_force = sec_area * PKII * (B_mat).transpose
    pub fn write_node_inner_force(&self, material_args: (Dtype, Dtype)) {
        let axial_force: Dtype = self.sec_area * self.PKII_stress(material_args)[0];
        let ax: Dtype = self.dx() / self.length_init();
        let ay: Dtype = self.dy() / self.length_init();
        let node_inner_force: [Dtype; 4] = [
            -axial_force * ax,
            -ay * axial_force,
            axial_force * ax,
            ay * axial_force,
        ];
        println!("inner_f: {:?}", &node_inner_force);
        for i in 0..2 {
            for j in 0..2 {
                self.nodes[i].forces[j] = node_inner_force[i * 2 as usize + j];
            }
        }
    }
    */
}

impl<'rod> K for Rod2D2NNL<'rod> {
    type Kmatrix = [[Dtype; 4]; 4];

    /// Cache stiffness matrix under linear analysis for rod element
    fn k(&mut self, material: (Dtype, Dtype)) -> &Self::Kmatrix
    where
        Self::Kmatrix: std::ops::Index<usize>,
    {
        if self.tangent_kmatrix.is_none() {
            self.tangent_kmatrix
                .get_or_insert(self.tangent_stiffness_mat(material))
        } else {
            self.tangent_kmatrix.as_ref().unwrap()
        }
    }

    fn k_printer(&self, n_exp: Dtype) {
        if self.tangent_kmatrix.is_none() {
            panic!(
                "!!! Rod2D2NNL#{}'s k mat is empty! call k() to calc it.",
                self.id
            );
        }

        print!("\nRod2D2NNL k{} =  (* 10^{})\n[", self.id, n_exp);
        for row in 0..4 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..4 {
                print!(
                    " {:>-10.6} ",
                    self.tangent_kmatrix.unwrap()[row][col]
                        / (10.0_f64.powf(n_exp as f64)) as Dtype
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

    /// Return rod elem's stiffness matrix's format string
    fn k_string(&self, n_exp: Dtype) -> String {
        let mut k_matrix = String::new();
        for row in 0..4 {
            if row == 0 {
                write!(k_matrix, "[[").expect("!!! Write rod k_mat to textfile failed!");
            } else {
                write!(k_matrix, " [").expect("!!! Write rod k_mat to textfile failed!");
            }
            for col in 0..4 {
                write!(
                    k_matrix,
                    " {:>-10.6} ",
                    self.tangent_kmatrix.unwrap()[row][col]
                        / (10.0_f64.powf(n_exp as f64)) as Dtype
                )
                .expect("!!! Write Rod2D2NNL kt_mat failed!");
            }
            if row == 1 {
                write!(k_matrix, "]]").expect("!!! Write rod k_mat to textfile failed!");
            } else {
                write!(k_matrix, "]\n").expect("!!! Write rod k_mat to textfile failed!");
            }
        }
        k_matrix
    }

    /// Get the stress at (x) inside the element
    fn strain(&self, _xyz: [Dtype; 3]) -> Vec<Dtype> {
        self.green_langrage_strain().to_vec()
    }

    /// Get the stress at (x) inside the element
    fn stress(&self, _xyz: [Dtype; 3], material: (Dtype, Dtype)) -> Vec<Dtype> {
        self.PKII_stress(material).to_vec()
    }

    /// Get element's info string
    fn info(&self) -> String {
        format!(
            "\n--------------------------------------------------------------------
Element_1D Info:\n\tId:     {}\n\tArea:   {}\n\tType:   Rod2D2NNL\n\tLen0:   {}\n\tLen:    {}\n\tNodes: {}\n\t       {}\n",
            self.id,
            self.sec_area,
            self.length_init(),
            self.length_current(),
            self.nodes[0],
            self.nodes[1]
        )
    }

    /// Get element's id number
    fn id(&self) -> usize {
        self.id
    }
}
