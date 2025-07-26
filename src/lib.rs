#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(test)]

extern crate nalgebra as na;
extern crate test;

mod calc;
mod dtty;
mod elem;
mod mesh;
mod node;
mod part;
mod port;
mod tool;

pub use calc::LinearEqs;
pub use dtty::{
    basic::{ADtype, Dtype, Jacobian2D, Jacobian3D},
    matrix::CompressedMatrix,
    sdata::Sdata,
};
pub use elem::{
    dim1::{beam::Beam1D2N, rod::Rod1D2N},
    dim2::{
        quadrila::Quad2D4N,
        rod::Rod2D2N,
        triangle::{Tri2D3N, Tri2D6N},
    },
};
pub use mesh::plane::Rectangle;
pub use node::{Node1D, Node2D, Node3D};
pub use part::{part1d::Part1D, part2d::Part2D};
pub use port::{Export, K};
pub use tool::*;

pub use std::collections::HashMap;

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

        assert_eq!([1.0 as Dtype], node1.coords);
        assert_eq!([1.0 as Dtype, 2.0 as Dtype], node2.coords);
        assert_eq!([1.0 as Dtype, 2.0 as Dtype, 3.0 as Dtype], node3.coords);
    }

    #[test]
    fn gen_elements() {
        let material: (Dtype, Dtype) = (1.0, 0.25);
        // 1D
        let node_a = Node1D::new(1, [0.0]);
        let node_b = Node1D::new(2, [1.0]);
        let node_c = Node1D::new(3, [3.0]);

        let rod1 = Rod1D2N::new(1, 1.0, [&node_a, &node_b], &material);
        let rod2 = Rod1D2N::new(2, 1.0, [&node_b, &node_c], &material);
        assert_eq!(1usize, rod1.id);
        assert_ne!(2usize, rod1.id);
        assert_eq!(2usize, rod2.id);

        // 2D
        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [0.0, 1.0]);
        let node3 = Node2D::new(3, [1.0, 0.0]);
        let node4 = Node2D::new(4, [1.0, 1.0]);
        let thick = 1.0;

        let tri1 = Tri2D3N::new(1, thick, [&node1, &node2, &node3], &material);
        let tri2 = Tri2D3N::new(2, thick, [&node4, &node2, &node3], &material);

        let rec1 = Quad2D4N::new(3, thick, [&node1, &node2, &node3, &node4], &material);

        assert_eq!(1usize, tri1.id);
        assert_ne!(2usize, tri1.id);
        assert_eq!(3usize, rec1.id);

        assert_ne!([10.0 as Dtype, 11.0 as Dtype], tri2.nodes[0].coords);
        assert_eq!([1.0 as Dtype, 1.0 as Dtype], tri2.nodes[0].coords);

        assert_eq!(vec![0.0, 0.0, 1.0], tri1.get_nodes_xcoords());
        assert_eq!(vec![0.0, 0.0, 1.0, 1.0], rec1.get_nodes_xcoords());
        assert_ne!(vec![0.0, 0.0, 1.0], tri1.get_nodes_ycoords());
        assert_ne!(vec![0.0, 0.0, 1.0, 1.0], rec1.get_nodes_ycoords());

        assert_eq!(0.5 as Dtype, tri1.area());
        assert_eq!(0.5 as Dtype, tri2.area());
        assert_eq!(1.0 as Dtype, rec1.area());
    }

    #[test]
    fn gen_parts() {
        let material = (1.0 as Dtype, 0.25 as Dtype);
        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [1.0, 0.0]);
        let node3 = Node2D::new(3, [1.0, 1.0]);
        let node4 = Node2D::new(4, [0.0, 1.0]);
        let thick = 1.0;

        let nodes = vec![node1, node2, node3, node4];
        let cplds = vec![vec![0, 1, 3], vec![1, 2, 3]];
        let mut tris: Vec<Tri2D3N> = tri2d3n_vec(thick, &nodes, &cplds, &material);

        let p1: Part2D<Tri2D3N, 4, 2, 3> = Part2D::new(1, &nodes, &mut tris, &cplds);
        assert_eq!(p1.elems[1].nodes[1].coords[1], 1.0);
        assert_ne!(p1.elems[1].nodes[1].coords[1], -1.0);

        let disp = vec![
            -1024., -1024., -1024., -1024., -1024., -1024., -1024., -1024.,
        ];
        let force = vec![0., 0., 0., 0., 0., 0., 0., 0.];
        let nodes_disp = p1.nodes_displacement();
        let nodes_force = p1.nodes_force();
        assert_ne!(disp, nodes_disp);
        assert_eq!(force, nodes_force);
    }

    #[test]
    fn mesh_rect_with_tri() {
        // set rect's width and height
        const W: Dtype = 1.0;
        const H: Dtype = 1.0;

        // number of nodes and freedom
        const R: usize = 2; //rows of nodes
        const C: usize = 2; //rows of nodes

        let rect_geo = Rectangle::new([0., 0.], [W, H]);
        let (coords, coupled_nodes) = rect_geo.mesh_with_tri2d3n(R, C);

        assert_eq!(coords, [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]);
        assert_eq!(coupled_nodes, [[0, 1, 2], [3, 2, 1]]);
        /*   2____3
         *   |\   |
         *   | \  |
         *   |  \ |
         *   |   \|
         *   0____1
         */
    }

    #[test]
    fn calc_tri_elem_k() {
        let material = (1.0 as Dtype, 0.25 as Dtype);

        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [1.0, 0.0]);
        let node3 = Node2D::new(3, [1.0, 1.0]);
        let node4 = Node2D::new(4, [0.0, 1.0]);
        let thick = 1.0;

        let mut tri1 = Tri2D3N::new(1, thick, [&node1, &node2, &node4], &material);
        let mut tri2 = Tri2D3N::new(2, thick, [&node3, &node4, &node2], &material);

        let k1: [[Dtype; 6]; 6] = tri1.k().recover();
        let k2: [[Dtype; 6]; 6] = tri2.k().recover();

        assert_eq!(k1, k2);
    }

    #[test]
    fn calc_quad_elem_k() {
        let material = (1.0 as Dtype, 0.25 as Dtype);

        let thick = 1.0;
        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [1.0, 0.0]);
        let node3 = Node2D::new(3, [1.0, 1.0]);
        let node4 = Node2D::new(4, [0.0, 1.0]);

        let mut quad1 = Quad2D4N::new(1, thick, [&node1, &node2, &node3, &node4], &material);

        let k1: [[Dtype; 8]; 8] = quad1.k().recover();
        assert_eq!(0.48888892 as Dtype, k1[0][0]);
    }

    #[bench]
    /// benchmark的结果是:393.49 ns +/- 36.23 ns/iter (Intel 8265U 插电) 20241205
    fn calc_quad_k_speed(b: &mut Bencher) {
        b.iter(|| calc_quad_elem_k());
    }

    #[bench]
    /// benchmark的结果是:229.13 ns +/- 06.28 ns/iter (Intel 8265U 插电) 20241205
    fn calc_tri_k_speed(b: &mut Bencher) {
        b.iter(|| calc_tri_elem_k());
    }
}
