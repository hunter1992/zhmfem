extern crate nalgebra as na;

use std::collections::HashMap;
use zhmfem::*;

fn main() {
    run();
}

fn run() {
    // set material parameters
    let thick = 1.0f64;
    let material = (1.0f64, 0.25f64);

    // input the node coordinates
    let coords: Vec<Vec<f64>> = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![1.0, 1.0],
        vec![0.0, 1.0],
    ];
    let zero_disp: Vec<usize> = vec![0, 1, 6];
    let force_index: Vec<usize> = vec![2, 4];
    let force_value: Vec<f64> = vec![-1.0, 1.0];
    let force_data: HashMap<usize, f64> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // transform points into nodes
    let nodes = nodes2d_vec(&coords, &zero_disp, force_data);

    // list nodes ids in one element
    let cpld: Vec<Vec<usize>> = vec![vec![1, 2, 4], vec![3, 4, 2]];

    // construct element by coupled nodes
    let mut tris: Vec<Tri2D3N> = tri2d3n_vec(thick, &nodes, &cpld);
    for i in tris.iter_mut() {
        println!("{}", i);
        i.k_printer(material);
    }

    // assemble global stiffness matrix
    let mut p1: Part2D<Tri2D3N, 4, 2, 3> = Part2D::new(1, tris, cpld);

    let mut solver = Solver::new(p1.disps(&nodes), p1.forces(&nodes));
    solver.static_kmat = Some(*p1.k(material));

    print_2darr("K", p1.k(material));
    print_1darr("qe", solver.get_disps());
    print_1darr("fe", solver.get_forces());
}
