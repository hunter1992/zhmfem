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

pub use calc::{
    panua::{pardiso, pardisoinit},
    solver::LinearEqs,
};
pub use dtty::{
    basic::{ADtype, Dtype, Jacobian2D, Jacobian3D},
    matrix::{CompressedMatrixCSR, CompressedMatrixSKS},
    sdata::NodeSData2D,
};
pub use elem::{
    dim1::rod::Rod1D2N,
    dim2::{quadrila::Quad2D4N, rod::Rod2D2N, triangle::Tri2D3N},
};
pub use mesh::plane::Rectangle;
pub use node::{Node1D, Node2D, Node3D};
pub use part::{part1d::Part1D, part2d::Part2D};
pub use port::{Export, StaticStiffness};
pub use tool::*;

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
        let material: [Dtype; 2] = [1.0, 0.25];

        // 1D
        let node_a = Node1D::new(1, [0.0]);
        let node_b = Node1D::new(2, [1.0]);
        let node_c = Node1D::new(3, [3.0]);

        let rod1 = Rod1D2N::new(1, 1.0, [&node_a, &node_b], material);
        let rod2 = Rod1D2N::new(2, 1.0, [&node_b, &node_c], material);
        assert_eq!(1usize, rod1.id);
        assert_ne!(2usize, rod1.id);
        assert_eq!(2usize, rod2.id);

        // 2D
        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [0.0, 1.0]);
        let node3 = Node2D::new(3, [1.0, 0.0]);
        let node4 = Node2D::new(4, [1.0, 1.0]);
        let thick = 1.0;

        let tri1 = Tri2D3N::new(1, thick, [&node1, &node2, &node3], material);
        let tri2 = Tri2D3N::new(2, thick, [&node4, &node2, &node3], material);

        // let rec1 = Quad2D4N::new(3, thick, [&node1, &node2, &node3, &node4], material);

        assert_eq!(1usize, tri1.id);
        assert_ne!(2usize, tri1.id);
        // assert_eq!(3usize, rec1.id);

        assert_ne!([10.0 as Dtype, 11.0 as Dtype], tri2.nodes[0].coords);
        assert_eq!([1.0 as Dtype, 1.0 as Dtype], tri2.nodes[0].coords);

        assert_eq!(vec![0.0, 0.0, 1.0], tri1.get_nodes_xcoord());
        // assert_eq!(vec![0.0, 0.0, 1.0, 1.0], rec1.get_nodes_xcoord());
        assert_ne!(vec![0.0, 0.0, 1.0], tri1.get_nodes_ycoord());
        // assert_ne!(vec![0.0, 0.0, 1.0, 1.0], rec1.get_nodes_ycoord());

        assert_eq!(0.5 as Dtype, tri1.area());
        assert_eq!(0.5 as Dtype, tri2.area());
        // assert_eq!(1.0 as Dtype, rec1.area());
    }

    #[test]
    fn gen_parts() {
        let ee = 1.0 as Dtype;
        let nu = 0.25 as Dtype;
        let thick = 1.0 as Dtype;
        let coords: Vec<[Dtype; 2]> = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];

        let cplds = vec![vec![0, 1, 3], vec![1, 2, 3]];

        let force_index: Vec<usize> = vec![2, 4];
        let force_value: Vec<Dtype> = vec![-1., 1.];
        let nodes = nodes2d_vec(&coords, &force_index, &force_value);

        let mut tris: Vec<Tri2D3N> = tri2d3n_vec(thick, &nodes, &cplds, [ee, nu]);

        let p1: Part2D<Tri2D3N, 4, 2, 3> = Part2D::new(1, &nodes, &mut tris, &cplds);
        assert_eq!(p1.elems[1].nodes[1].coords[1], 1.0);
        assert_ne!(p1.elems[1].nodes[1].coords[1], -1.0);

        let disp = vec![
            -1024., -1024., -1024., -1024., -1024., -1024., -1024., -1024.,
        ];
        let force = vec![0., 0., -1., 0., 1., 0., 0., 0.];
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
        let thick = 1.0;
        let material: [Dtype; 2] = [1.0, 0.25];

        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [1.0, 0.0]);
        let node3 = Node2D::new(3, [0.0, 1.0]);

        let nodes = [&node1, &node2, &node3];
        let tri1 = Tri2D3N::new(1, thick, nodes, material);
        let k: [[Dtype; 6]; 6] = tri1.calc_k();

        assert_eq!(0.7333333333333334, k[0][0]);
    }

    #[test]
    fn calc_quad_elem_k() {
        let thick = 1.0;
        let material: [Dtype; 2] = [1.0, 0.25];

        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [1.0, 0.0]);
        let node3 = Node2D::new(3, [1.0, 1.0]);
        let node4 = Node2D::new(4, [0.0, 1.0]);

        let nodes = [&node1, &node2, &node3, &node4];
        let quad1 = Quad2D4N::new(1, thick, nodes, material);
        let k1: [[Dtype; 8]; 8] = quad1.calc_k();

        assert_eq!(0.4888888888888888 as Dtype, k1[0][0]);
    }

    #[bench]
    /// benchmark的结果是:393.49 ns/iter (+/- 36.23) (Intel 8265U 插电) 20241205
    /// benchmark的结果是:285.15 ns/iter (+/- 32.51) (Intel 8265U 插电) 20251020
    fn calc_quad_k_speed(b: &mut Bencher) {
        b.iter(|| calc_quad_elem_k());
    }

    #[bench]
    /// benchmark的结果是:229.13 ns/iter (+/- 06.28) (Intel 8265U 插电) 20241205
    /// benchmark的结果是: 40.02 ns/iter (+/-  0.67) (Intel 8265U 插电) 20251020
    fn calc_tri_k_speed(b: &mut Bencher) {
        b.iter(|| calc_tri_elem_k());
    }
}
