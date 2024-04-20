use std::collections::HashMap;
use std::time::Instant;

use zhmfem::*;

fn main() {
    // set timing start
    let time_start = Instant::now();

    let material = (1.0 as Dtype, 0.25 as Dtype);

    let node1 = Node2D::new(1, [0., 0.]);
    let node2 = Node2D::new(2, [1., 0.]);
    let node3 = Node2D::new(3, [1., 1.]);
    let node4 = Node2D::new(4, [0., 1.]);

    // construct tri2d6n elements by coupled nodes
    let mut tri_vec: Vec<Tri2D6N> = vec![
        Tri2D6N::new(0, 1.0, [&node1, &node2, &node4]),
        Tri2D6N::new(1, 1.0, [&node3, &node4, &node2]),
    ];

    println!(
        "x: {}, {}, {}",
        tri_vec[0].xs()[0],
        tri_vec[0].xs()[1],
        tri_vec[0].xs()[2]
    );
    println!(
        "y: {}, {}, {}",
        tri_vec[0].ys()[0],
        tri_vec[0].ys()[1],
        tri_vec[0].ys()[2]
    );

    let b0 = tri_vec[1].geometry_mat_i(0, (1., 0.));
    print_2darr("B0", &b0, 0.0);
    let b1 = tri_vec[1].geometry_mat_i(1, (1., 0.));
    print_2darr("B1", &b1, 0.0);
    let b2 = tri_vec[1].geometry_mat_i(2, (1., 0.));
    print_2darr("B2", &b2, 0.0);
    let b3 = tri_vec[1].geometry_mat_i(3, (1., 0.));
    print_2darr("B3", &b3, 0.0);
    let b4 = tri_vec[1].geometry_mat_i(4, (1., 0.));
    print_2darr("B4", &b4, 0.0);
    let b5 = tri_vec[1].geometry_mat_i(5, (1., 0.));
    print_2darr("B5", &b5, 0.0);

    tri_vec[1].calc_k(material);
}
