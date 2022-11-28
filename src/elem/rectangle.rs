use super::triangle::Tri2D3N;
use crate::node::*;

pub struct Rec2D4N<'rect> {
    pub id: usize,
    pub thick: f64,
    pub nodes: [&'rect Node2D; 4],
}

impl Rec2D4N<'_> {
    /// generate a new Rec2D4N element
    pub fn new(id: usize, thick: f64, nodes: [&Node2D; 4]) -> Rec2D4N {
        Rec2D4N { id, thick, nodes }
    }

    /// get the x coords of nodes in Rec2D4N element
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
            "\nElement_2D Info:\n\tId:    {}\n\tArea:  {}\n\tType:  Rec2D4N
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
        let tri1: Tri2D3N = Tri2D3N {
            id: 0,
            thick: self.thick,
            nodes: [self.nodes[0], self.nodes[1], self.nodes[2]],
            k_matrix: None,
        };
        let tri2: Tri2D3N = Tri2D3N {
            id: 0,
            thick: self.thick,
            nodes: [self.nodes[3], self.nodes[1], self.nodes[2]],
            k_matrix: None,
        };
        let rect_area = tri1.area() + tri2.area();
        rect_area
    }
}
