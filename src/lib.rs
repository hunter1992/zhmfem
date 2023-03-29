#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(test)]

extern crate nalgebra as na;
extern crate test;

mod calc;
mod elem;
mod mesh;
mod node;
mod part;

pub use calc::Solver;
pub use elem::dim1::rod;
pub use elem::dim1::rod::Rod1D2N;
pub use elem::{dim2::quadrila::Quad2D4N, dim2::triangle::Tri2D3N};
pub use mesh::plane;
pub use na::*;
pub use node::*;
pub use part::Part2D;

use std::collections::HashMap;

type Jacobian2D = SMatrix<f64, 2, 2>;

pub trait K {
    type Kmatrix;
    fn k(&mut self, material: (f64, f64)) -> &Self::Kmatrix
    where
        Self::Kmatrix: std::ops::Index<usize>;
    fn k_printer(&self, n_exp: f64);
    fn k_string(&self, n_exp: f64) -> String;
}

pub fn print_1dvec<T>(name: &str, vec: &[T])
where
    T: std::fmt::Display,
{
    println!("{} =", name);
    print!("[[");
    for ele in vec.iter() {
        print!(" {:-7.4} ", &ele);
    }
    println!("]]\n");
}

pub fn print_2dvec<T>(name: &str, mat: &[Vec<T>])
where
    T: std::fmt::Display,
{
    println!("{} =", name);
    for row in 0..mat.len() {
        if row == 0 {
            print!("[[");
        } else {
            print!(" [");
        }
        for col in 0..mat[0].len() {
            print!(" {:-7.4} ", mat[row][col]);
        }
        if row == mat.len() - 1 {
            println!("]]\n");
        } else {
            println!("]");
        }
    }
}

pub fn print_1darr<T, const C: usize>(name: &str, arr: &[T; C])
where
    T: std::fmt::Display,
{
    println!("{} =", name);
    print!("[[");
    for c in 0..C {
        print!(" {:-7.4} ", arr[c]);
    }
    println!("]]\n");
}

pub fn print_2darr<T, const R: usize, const C: usize>(name: &str, arr: &[[T; C]; R])
where
    T: std::fmt::Display,
{
    println!("{} =", name);
    for r in 0..R {
        if r == 0 {
            print!("[[");
        } else {
            print!(" [");
        }
        for c in 0..C {
            print!(" {:-13.4} ", arr[r][c]);
        }
        if r == arr.len() - 1 {
            println!("]]\n");
        } else {
            println!("]");
        }
    }
}

pub fn nodes1d_vec(
    points: &[Vec<f64>],
    idx_0_disp: &[usize],
    force: &HashMap<usize, f64>,
) -> Vec<Node1D> {
    if points[0].len() != 1 {
        panic!(">>> Error from nodes1d_vec, the input points aren't 1D!");
    }

    let mut nodes: Vec<Node1D> = Vec::with_capacity(points.len());
    for (idx, coord) in points.iter().enumerate() {
        nodes.push(Node1D::new(idx + 1, [coord[0]]));
    }
    for idx in idx_0_disp.iter() {
        *nodes[idx / 1].disps[idx % 1].borrow_mut() = 0.0;
    }
    for (idx, &f) in force {
        *nodes[idx / 1].forces[idx % 1].borrow_mut() = f;
    }
    nodes
}

pub fn nodes2d_vec(
    points: &[Vec<f64>],
    idx_0_disp: &[usize],
    force: &HashMap<usize, f64>,
) -> Vec<Node2D> {
    if points[0].len() != 2 {
        panic!(">>> Error from nodes2d_vec, the input points aren't 2D!");
    }

    let mut nodes: Vec<Node2D> = Vec::with_capacity(points.len());
    for (idx, coord) in points.iter().enumerate() {
        nodes.push(Node2D::new(idx + 1, [coord[0], coord[1]]));
    }
    for idx in idx_0_disp.iter() {
        *nodes[idx / 2].disps[idx % 2].borrow_mut() = 0.0;
    }
    for (idx, &f) in force {
        *nodes[idx / 2].forces[idx % 2].borrow_mut() = f;
    }
    nodes
}

pub fn nodes3d_vec(points: &[Vec<f64>]) -> Vec<Node3D> {
    let mut nodes: Vec<Node3D> = Vec::new();
    for (idx, point) in points.iter().enumerate() {
        nodes.push(Node3D::new(idx + 1, [point[0], point[1], point[2]]));
    }
    nodes
}

/// Construct vector of rod1d2n elements
pub fn rod1d2n_vec<'rod>(
    sec_area: f64,
    nodes: &'rod [Node1D],
    couples: &[Vec<usize>],
) -> Vec<Rod1D2N<'rod>> {
    let mut rod1d2n: Vec<Rod1D2N> = Vec::new();
    for (ele_id, cpld) in couples.iter().enumerate() {
        rod1d2n.push(Rod1D2N::new(
            ele_id + 1,
            sec_area,
            [&nodes[cpld[0] - 1], &nodes[cpld[1] - 1]],
        ));
    }
    rod1d2n
}

/// Construct vector of tri2d3n elements
pub fn tri2d3n_vec<'tri>(
    thick: f64,
    nodes: &'tri [Node2D],
    couples: &[Vec<usize>],
) -> Vec<Tri2D3N<'tri>> {
    let mut tri2d3n: Vec<Tri2D3N> = Vec::new();
    for (ele_id, cpld) in couples.iter().enumerate() {
        tri2d3n.push(Tri2D3N::new(
            ele_id + 1,
            thick,
            [
                &nodes[cpld[0] - 1],
                &nodes[cpld[1] - 1],
                &nodes[cpld[2] - 1],
            ],
        ));
    }
    tri2d3n
}

/// Construct vector of quad2d4n elements
pub fn quad2d4n_vec<'rect>(
    thick: f64,
    nodes: &'rect [Node2D],
    couples: &[Vec<usize>],
) -> Vec<Quad2D4N<'rect>> {
    let mut rec2d4n: Vec<Quad2D4N> = Vec::new();
    for (ele_id, cpld) in couples.iter().enumerate() {
        rec2d4n.push(Quad2D4N::new(
            ele_id + 1,
            thick,
            [
                &nodes[cpld[0] - 1],
                &nodes[cpld[1] - 1],
                &nodes[cpld[2] - 1],
                &nodes[cpld[3] - 1],
            ],
        ))
    }
    rec2d4n
}

pub fn nonzero_index<'a, T: IntoIterator<Item = &'a f64>>(container: T) -> Vec<usize> {
    let idx: Vec<usize> = container
        .into_iter()
        .enumerate()
        .filter(|(_, &ele)| ele != 0.0)
        .map(|(idx, _)| idx)
        .collect();
    idx
}

#[cfg(test)]
mod testing {
    use super::*;
    use test::Bencher;

    #[test]
    fn gen_nodes() {
        let node1 = Node1D::new(1, [1.0]);
        let node2 = Node2D::new(2, [1.0, 2.0]);
        let node3 = Node3D::new(3, [1.0, 2.0, 3.0]);

        assert_eq!(1usize, node1.id);
        assert_ne!(3usize, node2.id);
        assert_eq!(3usize, node3.id);

        assert_eq!([1.0f64], node1.coord);
        assert_eq!([1.0f64, 2.0f64], node2.coord);
        assert_eq!([1.0f64, 2.0f64, 3.0f64], node3.coord);
    }

    #[test]
    fn gen_elements() {
        // 1D
        let node_a = Node1D::new(1, [0.0]);
        let node_b = Node1D::new(2, [1.0]);
        let node_c = Node1D::new(3, [3.0]);

        let rod1 = Rod1D2N::new(1, 1.0, [&node_a, &node_b]);
        let rod2 = Rod1D2N::new(2, 1.0, [&node_b, &node_c]);

        // 2D
        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [0.0, 1.0]);
        let node3 = Node2D::new(3, [1.0, 0.0]);
        let node4 = Node2D::new(4, [1.0, 1.0]);
        let thick = 1.0;

        let tri1 = Tri2D3N::new(1, thick, [&node1, &node2, &node3]);
        let tri2 = Tri2D3N::new(2, thick, [&node4, &node2, &node3]);

        let rec1 = Quad2D4N::new(3, thick, [&node1, &node2, &node3, &node4]);

        assert_eq!(1usize, tri1.id);
        assert_ne!(2usize, tri1.id);
        assert_eq!(3usize, rec1.id);

        assert_ne!([10.0f64, 11.0f64], tri2.nodes[0].coord);
        assert_eq!([1.0f64, 1.0f64], tri2.nodes[0].coord);

        assert_eq!(vec![0.0, 0.0, 1.0], tri1.xs());
        assert_eq!(vec![0.0, 0.0, 1.0, 1.0], rec1.xs());
        assert_ne!(vec![0.0, 0.0, 1.0], tri1.ys());
        assert_ne!(vec![0.0, 0.0, 1.0, 1.0], rec1.ys());

        assert_eq!(0.5f64, tri1.area());
        assert_eq!(0.5f64, tri2.area());
        assert_eq!(1.0f64, rec1.area());
    }

    #[test]
    fn gen_parts() {
        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [1.0, 0.0]);
        let node3 = Node2D::new(3, [1.0, -1.0]);
        let node4 = Node2D::new(4, [0.0, 1.0]);
        let thick = 1.0;

        let nodes = vec![node1, node2, node3, node4];
        let cplds = vec![vec![1, 2, 4], vec![2, 3, 4]];
        let mut tris: Vec<Tri2D3N> = tri2d3n_vec(thick, &nodes, &cplds);

        let p1: Part2D<Tri2D3N, 4, 2, 3> = Part2D::new(1, &nodes, &mut tris, &cplds);
        assert_eq!(p1.elems[1].nodes[1].coord[1], -1.0);
        assert_ne!(p1.elems[1].nodes[1].coord[1], 1.0);

        let disp = vec![
            -1024., -1024., -1024., -1024., -1024., -1024., -1024., -1024.,
        ];
        let force = vec![0., 0., 0., 0., 0., 0., 0., 0.];
        let nodes_disp = p1.disps();
        let nodes_force = p1.forces();
        assert_eq!(disp, nodes_disp);
        assert_eq!(force, nodes_force);
    }

    #[test]
    fn calc_elem_k() {
        let material = (1.0f64, 0.25f64);

        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [1.0, 0.0]);
        let node3 = Node2D::new(3, [1.0, 1.0]);
        let node4 = Node2D::new(4, [0.0, 1.0]);
        let thick = 1.0;

        let mut tri1 = Tri2D3N::new(1, thick, [&node1, &node2, &node4]);
        let mut tri2 = Tri2D3N::new(2, thick, [&node3, &node4, &node2]);

        let k1 = tri1.k(material);
        let k2 = tri2.k(material);

        assert_eq!(k1, k2);
    }

    #[test]
    fn mesh_rect_with_tri() {
        // set rect's width and height
        const W: f64 = 1.0;
        const H: f64 = 1.0;

        // number of nodes and freedom
        const R: usize = 2; //rows of nodes
        const C: usize = 2; //rows of nodes

        let rect_geo = plane::Rectangle::new([0., 0.], [W, H]);
        let (coords, coupled_nodes) = rect_geo.mesh_with_tri(R, C);

        assert_eq!(coords, [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]);
        assert_eq!(coupled_nodes, [[1, 2, 3], [4, 3, 2]]);
        /*   3____4
         *   |\   |
         *   | \  |
         *   |  \ |
         *   |___\|
         *   1    2
         */
    }

    #[bench]
    /// benchmark的结果是:277 +/- 15 ns/iter (Intel 8265U 插电)
    fn calc_elem_k_speed(b: &mut Bencher) {
        b.iter(|| calc_elem_k());
    }
}
