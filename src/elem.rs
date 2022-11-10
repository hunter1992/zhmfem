extern crate nalgebra as na;

use crate::node::*;
use na::*;
use std::fmt;

type Jacobian2x2f = SMatrix<f64, 2, 2>;

pub trait K {
    type Kmatrix;
    fn k(&mut self, material: (f64, f64, f64)) -> &Self::Kmatrix;
    fn k_printer(&mut self, material: (f64, f64, f64));
}

pub struct Triangle<'tri> {
    pub id: usize,
    pub nodes: [&'tri Node2D; 3],
    k_matrix: Option<[[f64; 6]; 6]>,
}

impl<'tri> Triangle<'tri> {
    /// generate a 2D triangle element
    pub fn new(id: usize, nodes: [&Node2D; 3]) -> Triangle {
        Triangle {
            id,
            nodes,
            k_matrix: None,
        }
    }

    /// get the x coords of nodes in tri element
    pub fn xs(&self) -> [f64; 3] {
        let mut x_list = [0.0; 3];
        for i in 0..3 {
            x_list[i] = self.nodes[i].coord[0];
        }
        x_list
    }

    /// get the y coords of nodes in tri element
    pub fn ys(&self) -> [f64; 3] {
        let mut y_list = [0.0; 3];
        for i in 0..3 {
            y_list[i] = self.nodes[i].coord[1];
        }
        y_list
    }

    /// calculate element stiffness matrix K
    /// return a 6x6 matrix, elements are f64
    fn calc_k(&self, material_args: (f64, f64, f64)) -> [[f64; 6]; 6] {
        let (ee, nu, t) = material_args;
        let elasticity_mat = SMatrix::<f64, 3, 3>::from([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - nu)],
        ]) * (ee / (1.0 - nu * nu));

        let x: [f64; 3] = self.xs();
        let y: [f64; 3] = self.ys();
        let dx21 = x[1] - x[0];
        let dx31 = x[2] - x[0];
        let dy21 = y[1] - y[0];
        let dy31 = y[2] - y[0];
        let jacobian = Jacobian2x2f::from([[dx21, dx31], [dy21, dy31]]);
        let det_j = jacobian.determinant();

        let h_mat = SMatrix::<f64, 3, 4>::from([
            [dy31, 0.0, -dx31],
            [-dy21, 0.0, dx21],
            [0.0, -dx31, dy31],
            [0.0, dx21, -dy21],
        ]) / (det_j.abs());

        let q_mat = SMatrix::<f64, 4, 6>::from([
            [-1., -1., 0., 0.],
            [0., 0., -1., -1.],
            [1., 0., 0., 0.],
            [0., 0., 1., 0.],
            [0., 1., 0., 0.],
            [0., 0., 0., 1.],
        ]);

        let b_mat = h_mat * q_mat;
        // Gauss integration, area of standard tri is 0.5
        let core = b_mat.transpose() * elasticity_mat * b_mat * det_j;
        let stiffness_matrix: [[f64; 6]; 6] = (0.5 * t * core).into();
        stiffness_matrix
    }

    pub fn area(&self) -> f64 {
        let x = self.xs();
        let y = self.ys();
        let tri_area = 0.5 * ((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])).abs();
        tri_area
    }
}

impl<'tri> K for Triangle<'tri> {
    type Kmatrix = [[f64; 6]; 6];

    /// cache stiffness matrix for element
    fn k(&mut self, material_args: (f64, f64, f64)) -> &Self::Kmatrix {
        if self.k_matrix.is_none() {
            self.k_matrix.get_or_insert(self.calc_k(material_args))
        } else {
            self.k_matrix.as_ref().unwrap()
        }
    }

    /// print element's stiffness matrix
    fn k_printer(&mut self, material_args: (f64, f64, f64)) {
        if self.k_matrix.is_none() {
            self.k_matrix = Some(self.calc_k(material_args));
        }
        print!("k{} = \n[", self.id);
        for row in 0..6 {
            if row == 0 {
                print!("[");
            } else {
                print!(" [")
            }
            for col in 0..6 {
                print!(" {:-9.6} ", self.k_matrix.unwrap()[row][col]);
            }
            if row == 5 {
                println!("]]");
            } else {
                println!("]")
            }
        }
        println!("");
    }
}

impl fmt::Display for Triangle<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\nElement_2D Info:\n\tId:     {}\n\tArea:   {}\n\tType:   Triangle
\tNodes: {}\n\t       {}\n\t       {}",
            self.id,
            self.area(),
            self.nodes[0],
            self.nodes[1],
            self.nodes[2]
        )
    }
}

pub struct Rectangle<'rect> {
    pub id: usize,
    pub nodes: [&'rect Node2D; 4],
}

impl Rectangle<'_> {
    /// generate a new rectangle element
    pub fn new(id: usize, nodes: [&Node2D; 4]) -> Rectangle {
        Rectangle { id, nodes }
    }

    /// get the x coords of nodes in rectangle element
    pub fn get_xs(&self) -> [f64; 4] {
        let mut xs = [0.0; 4];
        for i in 0..4 {
            xs[i] = self.nodes[i].coord[0];
        }
        xs
    }

    /// get the y coords of nodes in tri element
    pub fn get_ys(&self) -> [f64; 4] {
        let mut ys = [0.0; 4];
        for i in 0..4 {
            ys[i] = self.nodes[i].coord[1];
        }
        ys
    }

    pub fn info(&self) {
        println!(
            "\nElement_2D Info:\n\tId:    {}\n\tArea:  {}\n\tType:  Rectangle
\tNodes: {}\n\t       {}\n\t       {}\n\t       {}",
            self.id,
            self.area(),
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
            self.nodes[3]
        );
    }

    pub fn area(&self) -> f64 {
        let tri1: Triangle = Triangle {
            id: 0,
            nodes: [self.nodes[0], self.nodes[1], self.nodes[2]],
            k_matrix: None,
        };
        let tri2: Triangle = Triangle {
            id: 0,
            nodes: [self.nodes[3], self.nodes[1], self.nodes[2]],
            k_matrix: None,
        };
        let rect_area = tri1.area() + tri2.area();
        rect_area
    }
}
