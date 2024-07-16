use std::collections::HashMap;
use std::time::Instant;

use zhmfem::*;

fn main() {
    // set timing start
    let time_start = Instant::now();

    // set material parameters
    let thick = 1.0 as Dtype;
    let material = (1.0 as Dtype, 0.25 as Dtype); //Young's modulus & Poisson's ratio

    const W: Dtype = 1.0; // width
    const H: Dtype = 1.0; // height

    // number of nodes and freedom
    const R: usize = 2; // rows of nodes
    const C: usize = 2; // columns of nodes
    const M: usize = 3; // node num in single element
    const F: usize = 2; // freedom num in single node

    // construct the solid and mesh it
    let solid1 = plane::Rectangle::new([0.0 as Dtype, 0.0 as Dtype], [W, H]);
    let (coords, cpld) = solid1.mesh_with_tri2d6n(R, C);
    print_2dvec("coords", &coords, 0.0);
    print!("cpld=\n{:?}", cpld);

    //    let mut tri2d6n_vec:Vec<Tri2D6N> = tri2d6n_vec(thick, )

    /*
    let node1 = Node2D::new(1, [0.0, 0.0]);
    let node2 = Node2D::new(2, [1.0, 0.0]);
    let node3 = Node2D::new(3, [1.0, 1.0]);
    let node4 = Node2D::new(4, [0.0, 1.0]);
    let node5 = Node2D::new(5, [0.5, 0.0]);
    let node6 = Node2D::new(6, [1.0, 0.5]);
    let node7 = Node2D::new(7, [0.5, 1.0]);
    let node8 = Node2D::new(8, [0.0, 0.5]);
    let node9 = Node2D::new(9, [0.5, 0.5]);

    // construct tri2d6n elements by coupled nodes
    let mut tri_vec: Vec<Tri2D6N> = vec![
        Tri2D6N::new(0, 1.0, [&node1, &node2, &node3, &node5, &node6, &node9]),
        Tri2D6N::new(1, 1.0, [&node3, &node4, &node1, &node7, &node8, &node9]),
    ];

    let k_mat = tri_vec[1].calc_k(material);
    print_2darr("stiffness_matrix", &k_mat, 0.0);
    */
}
