use crate::{node::Node2D, Dtype, K};
use na::*;
use std::fmt::{self, Write};

pub struct Rod2D2NNL<'rod> {
    pub id: usize,
    pub sec_area: Dtype,
    pub nodes: [&'rod Node2D; 2],
    pub k_t_matrix: Option<[[Dtype; 4]; 4]>,
}

#[allow(non_snake_case)]
impl<'rod> Rod2D2NNL<'rod> {
    /// Generate a 2D Rod2D2N element
    pub fn new(id: usize, sec_area: Dtype, nodes: [&Node2D; 2]) -> Rod2D2NNL {
        Rod2D2NNL {
            id,
            sec_area,
            nodes,
            k_t_matrix: None,
        }
    }

    /// Get Node1's X coordinate under initial configuration
    pub fn X1(&self) -> Dtype {
        let X1: Dtype = self.nodes[0].coord[0];
        X1
    }

    /// Get Node1's x-coord displacement
    pub fn disp1_x(&self) -> Dtype {
        let disp1_x = *self.nodes[0].disps[0].borrow();
        disp1_x
    }

    /// Get Node1's x coordinate under current configuration
    pub fn x1(&self) -> Dtype {
        let x1: Dtype = self.X1() + self.disp1_x();
        x1
    }

    /// Get Node1's Y coordinate under initial configuration
    pub fn Y1(&self) -> Dtype {
        let Y1: Dtype = self.nodes[0].coord[1];
        Y1
    }

    /// Get Node1's y-coord displacement
    pub fn disp1_y(&self) -> Dtype {
        let disp1_y = *self.nodes[0].disps[1].borrow();
        disp1_y
    }

    /// Get Node1's y coordinate under current configuration
    pub fn y1(&self) -> Dtype {
        let y1: Dtype = self.Y1() + self.disp1_y();
        y1
    }

    /// Get Node2's X coordinate under initial configuration
    pub fn X2(&self) -> Dtype {
        let X2: Dtype = self.nodes[1].coord[0];
        X2
    }

    /// Get Node2's x-coord displacement
    pub fn disp2_x(&self) -> Dtype {
        let disp2_x = *self.nodes[1].disps[0].borrow();
        disp2_x
    }

    /// Get Node2's x coordinate under current configuration
    pub fn x2(&self) -> Dtype {
        let x2: Dtype = self.X2() + self.disp2_x();
        x2
    }

    /// Get Node2's Y coordinate under initial configuration
    pub fn Y2(&self) -> Dtype {
        let Y2: Dtype = self.nodes[1].coord[1];
        Y2
    }

    /// Get Node2's y-coord displacement
    pub fn disp2_y(&self) -> Dtype {
        let disp2_y = *self.nodes[1].disps[1].borrow();
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

    /// self displacement writer
    pub fn disp_writer(&self, disp_to_write: &[Dtype; 4]) {
        for i in 0..1 {
            for j in 0..1 {
                *self.nodes[i].disps[j].borrow_mut() = disp_to_write[i * 2 as usize + j];
            }
        }
    }

    /// Get strain matrix B of non-linear rod element
    /// B matrix  =  (1/L0) * [-ax, -ay, ax, ay], which is 1x4 mat
    /// and   ax  =  (x2 - x1) / L0,
    ///       ay  =  (y2 - y1) / L0
    pub fn strain_matrix(&self) -> [Dtype; 4] {
        let ax = (self.x2() - self.x1()) / self.length_init();
        let ay = (self.y2() - self.y1()) / self.length_init();
        [-ax, -ay, ax, ay]
    }

    /// Get Green-Lagrange strain of rod element
    pub fn strain_green(&self) -> Dtype {
        let l_init: Dtype = self.length_init();
        let l_current: Dtype = self.length_current();
        println!("  L0 = {}", l_init);
        println!("  L  = {}", l_current);
        let strain_green: Dtype =
            (l_current * l_current - l_init * l_init) / (2.0 as Dtype * l_init * l_init);
        strain_green
    }

    /// Get material stiffness matrix k_m of rod element
    pub fn k_m_stiffness_mat(&self, material_args: (Dtype, Dtype)) -> [[Dtype; 4]; 4] {
        let ee = material_args.0;
        let L0 = self.length_init();
        let b_mat = SMatrix::<Dtype, 1, 4>::from(self.strain_matrix());
        let km_mat: [[Dtype; 4]; 4] =
            (ee / L0 * self.sec_area * (b_mat.transpose() * b_mat)).into();
        km_mat
    }

    /// Get geometry stiffness matrix k_g of rod element
    pub fn k_g_stiffness_mat(&self, material_args: (Dtype, Dtype)) -> [[Dtype; 4]; 4] {
        let ee = material_args.0;
        let arg: Dtype = ee * self.sec_area * self.strain_green() / self.length_init();
        let base_mat = SMatrix::<Dtype, 4, 4>::from([
            [1.0, 0.0, -1.0, 0.0],
            [0., 1., 0., -1.],
            [-1., 0., 1., 0.],
            [0., -1., 0., 1.],
        ]);
        let kg_mat: [[Dtype; 4]; 4] = (arg * base_mat).into();
        kg_mat
    }
}

pub struct Rod2D2N<'rod> {
    pub id: usize,
    pub sec_area: Dtype,
    pub nodes: [&'rod Node2D; 2],
    pub k_matrix: Option<[[Dtype; 4]; 4]>,
}

impl<'rod> Rod2D2N<'rod> {
    /// Generate a 2D Rod2D2N element
    pub fn new(id: usize, sec_area: Dtype, nodes: [&Node2D; 2]) -> Rod2D2N {
        Rod2D2N {
            id,
            sec_area,
            nodes,
            k_matrix: None,
        }
    }

    /// Get id number
    pub fn get_id(&self) -> usize {
        let id_num: usize = self.id;
        id_num
    }

    /// Get the x-coords of nodes in rod element
    pub fn xs(&self) -> [Dtype; 2] {
        let mut x_list = [0.0; 2];
        for i in 0..2 {
            x_list[i] = self.nodes[i].coord[0];
        }
        x_list
    }

    /// Get the y-coords of nodes in rod element
    pub fn ys(&self) -> [Dtype; 2] {
        let mut y_list = [0.0; 2];
        for i in 0..2 {
            y_list[i] = self.nodes[i].coord[1];
        }
        y_list
    }

    /// Get rod element length
    pub fn length(&self) -> (Dtype, Dtype, Dtype) {
        let dx = self.xs()[1] - self.xs()[0];
        let dy = self.ys()[1] - self.ys()[0];
        let l: Dtype = (dx * dx + dy * dy).sqrt();
        let cos: Dtype = ((dx as f64) / (l as f64)) as Dtype;
        let sin: Dtype = ((dy as f64) / (l as f64)) as Dtype;
        (l, cos, sin)
    }

    /// Get nodes' disps vector
    pub fn disps(&self) -> [Dtype; 4] {
        let mut disps = [0.0; 4];
        for idx in 0..2 {
            disps[2 * idx] = *self.nodes[idx].disps[0].borrow();
            disps[2 * idx + 1] = *self.nodes[idx].disps[1].borrow();
        }
        disps
    }

    /// Get nodes' force vector
    pub fn forces(&self) -> [Dtype; 4] {
        let mut forces = [0.0; 4];
        for idx in 0..2 {
            forces[2 * idx] = *self.nodes[idx].forces[0].borrow();
            forces[2 * idx + 1] = *self.nodes[idx].forces[1].borrow();
        }
        forces
    }

    /// Get transformation matrix
    pub fn trans_mat(&self) -> [[Dtype; 2]; 4] {
        let (_, cos, sin) = self.length();
        [[cos, 0.0], [sin, 0.0], [0.0, cos], [0.0, sin]]
    }

    /// Calculate element stiffness matrix K
    /// Return a 4x4 matrix, elements are Dtype
    fn calc_k(&self, material_args: (Dtype, Dtype)) -> [[Dtype; 4]; 4] {
        println!(
            "\n>>> Calculating Rod2D2N(#{})'s stiffness matrix k{}",
            self.id, self.id
        );
        let (ee, _nu) = material_args;
        let trans_mat = SMatrix::<Dtype, 2, 4>::from(self.trans_mat());
        let local_stiffness_mat = SMatrix::<Dtype, 2, 2>::from([[1.0, -1.0], [-1.0, 1.0]])
            * (ee * self.sec_area / self.length().0);
        let global_stiffness_mat: [[Dtype; 4]; 4] =
            (trans_mat.transpose() * local_stiffness_mat * trans_mat).into();
        global_stiffness_mat
    }

    /// Get element's strain vector, a scale in rod element
    fn calc_strain(&self) -> [Dtype; 3] {
        let unit = 1.0 / self.length().0;
        let b_mat = SMatrix::<Dtype, 1, 2>::from([-unit, unit]);
        let trans_mat = SMatrix::<Dtype, 2, 4>::from(self.trans_mat());
        let nodes_disps = SMatrix::<Dtype, 4, 1>::from(self.disps());
        let strain: [Dtype; 1] = (b_mat * trans_mat * nodes_disps).into();
        [strain[0], 0.0, 0.0]
    }

    /// Get element's strain vector, a scale in rod element
    fn calc_stress(&self, material_args: (Dtype, Dtype)) -> [Dtype; 3] {
        let (ee, _nu) = material_args;
        let stress: [Dtype; 1] = [ee * self.calc_strain()[0]];
        [stress[0], 0.0, 0.0]
    }
}

impl<'rod> K for Rod2D2N<'rod> {
    type Kmatrix = [[Dtype; 4]; 4];

    /// Cache stiffness matrix for rod element
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

    /// Print rod element's stiffness matrix
    fn k_printer(&self, n_exp: Dtype) {
        if self.k_matrix.is_none() {
            panic!(
                "!!! Rod2D2N#{}'s k mat is empty! call k() to calc it.",
                self.id
            );
        }

        print!("\nRod2D2N k{} =  (* 10^{})\n[", self.id, n_exp);
        for row in 0..4 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..4 {
                print!(
                    " {:>-10.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
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
                write!(k_matrix, "[[").expect("!!! Write tri k_mat failed!");
            } else {
                write!(k_matrix, " [").expect("!!! Write tri k_mat failed!");
            }
            for col in 0..4 {
                write!(
                    k_matrix,
                    " {:>-10.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
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

    /// Get the stress at (x) inside the element
    fn strain(&self, _xyz: [Dtype; 3]) -> Vec<Dtype> {
        self.calc_strain().to_vec()
    }

    /// Get the stress at (x) inside the element
    fn stress(&self, _xyz: [Dtype; 3], material: (Dtype, Dtype)) -> Vec<Dtype> {
        self.calc_stress(material).to_vec()
    }

    /// Get element's info string
    fn info(&self) -> String {
        format!(
            "\n--------------------------------------------------------------------
Element_1D Info:\n\tId:     {}\n\tArea:   {}\n\tLong:   {}\n\tType:   Rod2D2N\n\tNodes: {}\n\t       {}\n",
            self.id,
            self.sec_area,
            self.length().0,
            self.nodes[0],
            self.nodes[1]
        )
    }

    /// Get element's id number
    fn id(&self) -> usize {
        self.id
    }
}

impl fmt::Display for Rod2D2N<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,"\n--------------------------------------------------------------------
Element_1D Info:\n\tId:     {}\n\tArea:   {}\n\tLong:   {}\n\tType:   Rod2D2N\n\tNodes: {}\n\t       {}\n",
self.id, self.sec_area, self.length().0, self.nodes[0], self.nodes[1])
    }
}
