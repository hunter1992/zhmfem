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

pub use calc::LinearEqs;
pub use elem::dim1::{beam::Beam1D2N, rod::Rod1D2N};
pub use elem::dim2::linear::{
    quadrila::Quad2D4N,
    rod::Rod2D2N,
    triangle::{Tri2D3N, Tri2D6N},
};
pub use elem::dim2::nonlinear::rod::Rod2D2NNL;
pub use mesh::plane;
pub use na::*;
pub use node::*;
pub use part::{part1d::Part1D, part2d::Part2D};

use std::collections::HashMap;

pub type Dtype = f32;
pub type Jacobian2D = SMatrix<Dtype, 2, 2>;

/// K trait generate element's stiffness matrix under linear analysis.
/// Output stress/strain vector at some point in element using Vector,
/// Kmatrix is the element's stiffness matrix to get.
pub trait K {
    type Kmatrix;

    fn k(&mut self, material: (Dtype, Dtype)) -> &Self::Kmatrix
    where
        Self::Kmatrix: std::ops::Index<usize>;
    fn k_printer(&self, n_exp: Dtype);
    fn k_string(&self, n_exp: Dtype) -> String;

    fn strain(&self, xyz: [Dtype; 3]) -> Vec<Dtype>;
    fn stress(&self, xyz: [Dtype; 3], material: (Dtype, Dtype)) -> Vec<Dtype>;

    fn info(&self) -> String;
    fn id(&self) -> usize;
}

pub trait Export {
    fn txt_writer(&self, target_file: &str) -> std::io::Result<bool>;
    fn vtk_writer(&self, target_file: &str) -> std::io::Result<bool>;
}

/// 用slave矩阵填充master矩阵中的一部分，填充起始位置由sp决定
#[inline]
pub fn matrix_block_fill<const R1: usize, const C1: usize, const R2: usize, const C2: usize>(
    master: &mut [[Dtype; C1]; R1],
    slave: &[[Dtype; C2]; R2],
    sp: (usize, usize),
) {
    for row in 0..slave.len() {
        for col in 0..slave[0].len() {
            master[sp.0 + row][sp.1 + col] = slave[row][col];
        }
    }
}

/// 用slave矩阵追加master矩阵中的一部分，追加起始位置由sp决定
#[inline]
pub fn matrix_block_append<const R1: usize, const C1: usize, const R2: usize, const C2: usize>(
    master: &mut [[Dtype; C1]; R1],
    slave: &[[Dtype; C2]; R2],
    sp: (usize, usize),
) {
    for row in 0..slave.len() {
        for col in 0..slave[0].len() {
            master[sp.0 + row][sp.1 + col] += slave[row][col];
        }
    }
}

#[inline]
pub fn print_1dvec(name: &str, vec: &[Dtype], n_exp: Dtype) {
    println!("\n{} =", name);
    print!("[[");
    for &ele in vec.iter() {
        print!(" {:-9.4} ", ele / (10.0_f64.powf(n_exp as f64)) as Dtype);
    }
    println!("]]\n");
}

#[inline]
pub fn print_2dvec(name: &str, mat: &[Vec<Dtype>], n_exp: Dtype) {
    println!("\n{} =", name);
    for row in 0..mat.len() {
        if row == 0 {
            print!("[[");
        } else {
            print!(" [");
        }
        for col in 0..mat[0].len() {
            print!(
                " {:-9.4} ",
                mat[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
            );
        }
        if row == mat.len() - 1 {
            println!("]]\n");
        } else {
            println!("]");
        }
    }
}

/// formatted print 2d array with scientific form
#[inline]
pub fn print_2darr<const R: usize, const C: usize>(
    name: &str,
    arr: &[[Dtype; C]; R],
    n_exp: Dtype,
) {
    println!("\n{} = (* 10^{})", name, n_exp);
    for r in 0..R {
        if r == 0 {
            print!("[[");
        } else {
            print!(" [");
        }
        for c in 0..C {
            print!(
                " {:-10.6} ",
                arr[r][c] / (10.0_f64.powf(n_exp as f64)) as Dtype
            );
        }
        if r == arr.len() - 1 {
            println!("]]\n");
        } else {
            println!("]");
        }
    }
}

pub fn nodes1d_vec(points: &[Vec<Dtype>], force: &HashMap<usize, Dtype>) -> Vec<Node1D> {
    if points[0].len() != 1_usize {
        panic!(">>> Error from nodes1d_vec, the input points aren't 1D!");
    }

    let mut nodes: Vec<Node1D> = Vec::with_capacity(points.len());
    for (idx, coord) in points.iter().enumerate() {
        nodes.push(Node1D::new(idx, [coord[0]]));
    }
    for (idx, &f) in force {
        *nodes[idx / 1].forces[idx % 1].borrow_mut() = f;
    }
    nodes
}

pub fn nodes2d_vec(
    points: &[Vec<Dtype>],
    force: &HashMap<usize, Dtype>,
    nonlinear_or_not: bool,
) -> Vec<Node2D> {
    if points[0].len() != 2 {
        panic!(">>> Error from nodes2d_vec, the input points aren't 2D!");
    }

    let mut nodes: Vec<Node2D> = Vec::with_capacity(points.len());
    for (idx, coord) in points.iter().enumerate() {
        nodes.push(Node2D::new(idx, [coord[0], coord[1]]));
    }

    if nonlinear_or_not { //如果是线性分析，不区分节点的内力外力；非线性需要区分
    } else {
        for (idx, &f) in force {
            *nodes[idx / 2].forces[idx % 2].borrow_mut() = f;
        }
    }
    nodes
}

pub fn nodes3d_vec(points: &[Vec<Dtype>]) -> Vec<Node3D> {
    let mut nodes: Vec<Node3D> = Vec::new();
    for (idx, point) in points.iter().enumerate() {
        nodes.push(Node3D::new(idx, [point[0], point[1], point[2]]));
    }
    nodes
}

/// Construct vector of rod1d2n elements
pub fn rod1d2n_vec<'rod>(
    sec_areas: &'rod [Dtype],
    nodes: &'rod [Node1D],
    coupled_nodes: &[Vec<usize>],
) -> Vec<Rod1D2N<'rod>> {
    let mut rod1d2n: Vec<Rod1D2N> = Vec::with_capacity(sec_areas.len());
    for (ele_id, nodes_id_pair) in coupled_nodes.iter().enumerate() {
        rod1d2n.push(Rod1D2N::new(
            ele_id,
            sec_areas[ele_id],
            [&nodes[nodes_id_pair[0]], &nodes[nodes_id_pair[1]]],
        ));
    }
    rod1d2n
}

/// Construct vector of rod2d2n elements
pub fn rod2d2n_vec<'rod>(
    sec_areas: &'rod [Dtype],
    nodes: &'rod [Node2D],
    coupled_nodes: &[Vec<usize>],
) -> Vec<Rod2D2N<'rod>> {
    let mut rod2d2n: Vec<Rod2D2N> = Vec::with_capacity(sec_areas.len());
    for (ele_id, nodes_id_pair) in coupled_nodes.iter().enumerate() {
        rod2d2n.push(Rod2D2N::new(
            ele_id,
            sec_areas[ele_id],
            [&nodes[nodes_id_pair[0]], &nodes[nodes_id_pair[1]]],
        ));
    }
    rod2d2n
}

/// Construct vector of rod2d2n_nonlinear elements
pub fn rod2d2n_nonlinear_vec<'rod>(
    sec_areas: &'rod [Dtype],
    nodes: &'rod [Node2D],
    coupled_nodes: &[Vec<usize>],
) -> Vec<Rod2D2NNL<'rod>> {
    let mut rod2d2nnl: Vec<Rod2D2NNL> = Vec::with_capacity(sec_areas.len());
    for (ele_id, nodes_id_pair) in coupled_nodes.iter().enumerate() {
        rod2d2nnl.push(Rod2D2NNL::new(
            ele_id,
            sec_areas[ele_id],
            [&nodes[nodes_id_pair[0]], &nodes[nodes_id_pair[1]]],
        ));
    }
    rod2d2nnl
}

/// Construct vector of beam1d2n elements
pub fn beam1d2n_vec<'beam>(
    moi: Dtype,
    sec_area: Dtype,
    nodes: &'beam [Node2D],
    coupled_nodes: &[Vec<usize>],
) -> Vec<Beam1D2N<'beam>> {
    let mut beam1d2n_vec: Vec<Beam1D2N> = Vec::new();
    for (ele_id, nodes_id_pair) in coupled_nodes.iter().enumerate() {
        beam1d2n_vec.push(Beam1D2N::new(
            ele_id,
            moi,
            sec_area,
            [&nodes[nodes_id_pair[0]], &nodes[nodes_id_pair[1]]],
        ));
    }
    beam1d2n_vec
}

/// Construct vector of tri2d3n elements
pub fn tri2d3n_vec<'tri>(
    thick: Dtype,
    nodes: &'tri [Node2D],
    couples: &[Vec<usize>],
) -> Vec<Tri2D3N<'tri>> {
    let mut tri2d3n: Vec<Tri2D3N> = Vec::new();
    for (ele_id, cpld) in couples.iter().enumerate() {
        tri2d3n.push(Tri2D3N::new(
            ele_id,
            thick,
            [&nodes[cpld[0]], &nodes[cpld[1]], &nodes[cpld[2]]],
        ));
    }
    tri2d3n
}

/// Construct vector of quad2d4n elements
pub fn quad2d4n_vec<'rect>(
    thick: Dtype,
    nodes: &'rect [Node2D],
    couples: &[Vec<usize>],
) -> Vec<Quad2D4N<'rect>> {
    let mut rec2d4n: Vec<Quad2D4N> = Vec::new();
    for (ele_id, cpld) in couples.iter().enumerate() {
        rec2d4n.push(Quad2D4N::new(
            ele_id,
            thick,
            [
                &nodes[cpld[0]],
                &nodes[cpld[1]],
                &nodes[cpld[2]],
                &nodes[cpld[3]],
            ],
        ))
    }
    rec2d4n
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

        assert_eq!([1.0 as Dtype], node1.coord);
        assert_eq!([1.0 as Dtype, 2.0 as Dtype], node2.coord);
        assert_eq!([1.0 as Dtype, 2.0 as Dtype, 3.0 as Dtype], node3.coord);
    }

    #[test]
    fn gen_elements() {
        // 1D
        let node_a = Node1D::new(1, [0.0]);
        let node_b = Node1D::new(2, [1.0]);
        let node_c = Node1D::new(3, [3.0]);

        let rod1 = Rod1D2N::new(1, 1.0, [&node_a, &node_b]);
        let rod2 = Rod1D2N::new(2, 1.0, [&node_b, &node_c]);
        assert_eq!(1usize, rod1.id);
        assert_ne!(2usize, rod1.id);
        assert_eq!(2usize, rod2.id);

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

        assert_ne!([10.0 as Dtype, 11.0 as Dtype], tri2.nodes[0].coord);
        assert_eq!([1.0 as Dtype, 1.0 as Dtype], tri2.nodes[0].coord);

        assert_eq!(vec![0.0, 0.0, 1.0], tri1.xs());
        assert_eq!(vec![0.0, 0.0, 1.0, 1.0], rec1.xs());
        assert_ne!(vec![0.0, 0.0, 1.0], tri1.ys());
        assert_ne!(vec![0.0, 0.0, 1.0, 1.0], rec1.ys());

        assert_eq!(0.5 as Dtype, tri1.area());
        assert_eq!(0.5 as Dtype, tri2.area());
        assert_eq!(1.0 as Dtype, rec1.area());
    }

    #[test]
    fn gen_parts() {
        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [1.0, 0.0]);
        let node3 = Node2D::new(3, [1.0, 1.0]);
        let node4 = Node2D::new(4, [0.0, 1.0]);
        let thick = 1.0;

        let nodes = vec![node1, node2, node3, node4];
        let cplds = vec![vec![0, 1, 3], vec![1, 2, 3]];
        let zero_disps_idx = vec![0, 1, 3, 6];
        let mut tris: Vec<Tri2D3N> = tri2d3n_vec(thick, &nodes, &cplds);

        let p1: Part2D<Tri2D3N, 4, 2, 3> =
            Part2D::new(1, &nodes, &mut tris, &cplds, &zero_disps_idx);
        assert_eq!(p1.elems[1].nodes[1].coord[1], 1.0);
        assert_ne!(p1.elems[1].nodes[1].coord[1], -1.0);

        let disp = vec![
            -1024., -1024., -1024., -1024., -1024., -1024., -1024., -1024.,
        ];
        let force = vec![0., 0., 0., 0., 0., 0., 0., 0.];
        let nodes_disp = p1.disps();
        let nodes_force = p1.forces();
        assert_ne!(disp, nodes_disp);
        assert_eq!(force, nodes_force);
    }

    #[test]
    fn calc_tri_elem_k() {
        let material = (1.0 as Dtype, 0.25 as Dtype);

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
    fn calc_quad_elem_k() {
        let material = (1.0 as Dtype, 0.25 as Dtype);

        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [1.0, 0.0]);
        let node3 = Node2D::new(3, [1.0, 1.0]);
        let node4 = Node2D::new(4, [0.0, 1.0]);
        let thick = 1.0;

        let mut quad1 = Quad2D4N::new(1, thick, [&node1, &node2, &node3, &node4]);

        let k1 = quad1.k(material);
        assert_eq!(0.48888892 as Dtype, k1[0][0]);
    }

    #[test]
    fn mesh_rect_with_tri() {
        // set rect's width and height
        const W: Dtype = 1.0;
        const H: Dtype = 1.0;

        // number of nodes and freedom
        const R: usize = 2; //rows of nodes
        const C: usize = 2; //rows of nodes

        let rect_geo = plane::Rectangle::new([0., 0.], [W, H]);
        let (coords, coupled_nodes) = rect_geo.mesh_with_tri(R, C);

        assert_eq!(coords, [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]);
        assert_eq!(coupled_nodes, [[0, 1, 2], [3, 2, 1]]);
        /*   2____3
         *   |\   |
         *   | \  |
         *   |  \ |
         *   |___\|
         *   0    1
         */
    }

    #[bench]
    /// benchmark的结果是:277 +/- 15 ns/iter (Intel 8265U 插电)
    fn calc_elem_k_speed(b: &mut Bencher) {
        b.iter(|| calc_quad_elem_k());
    }
}
