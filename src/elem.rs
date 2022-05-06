use crate::node::*;

pub struct Triangle<'tri> {
    pub id: usize,
    pub nodes: [&'tri Node2D; 3],
    k_matrix: Option<Vec<Vec<f64>>>,
}

impl <'tri> Triangle<'tri> {
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
    pub fn k(&mut self, args: (f64, f64, f64)) -> Vec<Vec<f64>> {
        match self.k_matrix {
            Some(k_mat) => &k_mat,
            None => {
                let k_mat = self.calc_k(args);
                self.k_matrix = Some(k_mat);
                &self.k_matrix.unwrap()
            },
        }
    }

    /// calculate element stiffness matrix K
    fn calc_k(&self, args: (f64, f64, f64)) -> Vec<Vec<f64>> {
        let (ee, nu, t) = args;
        let xs: [f64; 3] = self.xs();
        let ys: [f64; 3] = self.ys();

        let idx = |x| x % 3;
        let b = c![ys[idx(i+1) as usize] - ys[idx(i+2) as usize],
                   for i in 0..3];
        let c = c![xs[idx(i+2) as usize] - xs[idx(i+1) as usize],
                   for i in 0..3];

        let k = |i: usize, j: usize| {
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

        let pick_k = |num: usize| {
            if num % 2 == 0 {
                (num / 2) as usize
            } else {
                ((num - 1) / 2) as usize
            }
        };
        let pick_kij = |num: usize| (num % 2) as usize;
        let coef = ee * t / (4.0 * (1.0 - nu * nu) * self.area());
        let k_mat = c![c![ coef * k(pick_k(i),
                           pick_k(j))[pick_kij(i)][pick_kij(j)],
                           for i in 0..6 ],
                           for j in 0..6 ];
        k_mat
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
