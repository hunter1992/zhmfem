use crate::node::*;

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

    /// cache stiffness matrix for element
    pub fn k(&mut self, material_args: (f64, f64, f64)) -> &[[f64; 6]; 6] {
        if self.k_matrix.is_none() {
            self.k_matrix.get_or_insert(self.calc_k(material_args))
        } else {
            self.k_matrix.as_ref().unwrap()
        }
    }

    /// calculate element stiffness matrix K
    /// return a 6x6 matrix, elements are f64
    fn calc_k(&self, material_args: (f64, f64, f64)) -> [[f64; 6]; 6] {
        println!("--->Calculating triangle[{}]'s stiffness matrix...", self.id);
        let (ee, nu, t) = material_args;
        let xs: [f64; 3] = self.xs();
        let ys: [f64; 3] = self.ys();

        let idx = |x| x % 3;
        let b = [
            ys[idx(0 + 1) as usize] - ys[idx(0 + 2) as usize],
            ys[idx(1 + 1) as usize] - ys[idx(1 + 2) as usize],
            ys[idx(2 + 1) as usize] - ys[idx(2 + 2) as usize],
        ];
        let c = [
            xs[idx(0 + 2) as usize] - xs[idx(0 + 1) as usize],
            xs[idx(1 + 2) as usize] - xs[idx(1 + 1) as usize],
            xs[idx(2 + 2) as usize] - xs[idx(2 + 1) as usize],
        ];

        let k_mat = |i: usize, j: usize| {
            [
                [
                    b[i] * b[j] + 0.5 * (1.0 - nu) * c[i] * c[j],
                    nu * b[i] * c[j] + 0.5 * (1.0 - nu) * c[i] * b[j],
                ],
                [
                    nu * c[i] * b[j] + 0.5 * (1.0 - nu) * b[i] * c[j],
                    c[i] * c[j] + 0.5 * (1.0 - nu) * b[i] * b[j],
                ],
            ]
        };
        let pick_kij = |num: usize| {
            if num % 2 == 0 {
                (num / 2) as usize
            } else {
                ((num - 1) / 2) as usize
            }
        };
        let ij = |num: usize| (num % 2) as usize;
        let coef = ee * t / (4.0 * (1.0 - nu * nu) * self.area());
        let k = |i: usize, j: usize| coef * k_mat(pick_kij(i), pick_kij(j))[ij(i)][ij(j)];
        let stiffness_matrix = [
            [k(0, 0), k(0, 1), k(0, 2), k(0, 3), k(0, 4), k(0, 5)],
            [k(1, 0), k(1, 1), k(1, 2), k(1, 3), k(1, 4), k(1, 5)],
            [k(2, 0), k(2, 1), k(2, 2), k(2, 3), k(2, 4), k(2, 5)],
            [k(3, 0), k(3, 1), k(3, 2), k(3, 3), k(3, 4), k(3, 5)],
            [k(4, 0), k(4, 1), k(4, 2), k(4, 3), k(4, 4), k(4, 5)],
            [k(5, 0), k(5, 1), k(5, 2), k(5, 3), k(5, 4), k(5, 5)],
        ];
        stiffness_matrix
    }

    pub fn k_printer(&mut self, material_args: (f64, f64, f64)) {
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

    pub fn info(&self) {
        println!(
            "\nElement_2D Info:\n\tId:    {}\n\tArea:  {}\n\tType:  Triangle
\tCoord: {:?}\n\t       {:?}\n\t       {:?}",
            self.id,
            self.area(),
            self.nodes[0],
            self.nodes[1],
            self.nodes[2]
        );
        println!("  If you wanna:");
        println!("        see the stiffness matrix, use:");
        println!("               tri-elem_name.k_printer((ee, nu, t))");
        println!("        or, just get the stiffness matrix:");
        println!("               tri-elem_name.k((ee, nu, t))");
        print!("\n");
    }

    pub fn area(&self) -> f64 {
        let x1 = self.nodes[0].coord[0];
        let y1 = self.nodes[0].coord[1];
        let x2 = self.nodes[1].coord[0];
        let y2 = self.nodes[1].coord[1];
        let x3 = self.nodes[2].coord[0];
        let y3 = self.nodes[2].coord[1];
        let tri_area = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)).abs();
        tri_area
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
\tCoord: {:?}\n\t       {:?}\n\t       {:?}\n\t       {:?}",
            self.id,
            self.area(),
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
            self.nodes[3]
        );
    }

    pub fn area(&self) -> f64 {
        let tri1 = Triangle {
            id: 0,
            nodes: [self.nodes[0], self.nodes[1], self.nodes[2]],
            k_matrix: None,
        };
        let tri2 = Triangle {
            id: 0,
            nodes: [self.nodes[3], self.nodes[1], self.nodes[2]],
            k_matrix: None,
        };
        let rect_area = tri1.area() + tri2.area();
        rect_area
    }
}
