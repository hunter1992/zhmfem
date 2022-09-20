//extern crate nalgebra as na;

use zhmfem::*;

fn main() {
    run();
}

fn run() {
    // set material parameters
    let ee = 1.0f64;
    let nu = 0.25f64;
    let t = 1.0f64;
    let material = (ee, nu, t);

    // input the points coordinates
    let points: Vec<Vec<f64>> = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![1.0, 1.0],
        vec![0.0, 1.0],
    ];

    // transform points into nodes
    let nodes = nodes2d_vec(&points);

    // list nodes ids in one element
    let coupled_nodes: Vec<Vec<usize>> = vec![vec![1, 2, 4], vec![2, 3, 4]];

    // construct element by coupled nodes
    let mut tris: Vec<Triangle> = tri2d3n_vec(&nodes, &coupled_nodes);
    println!("{}", tris[0]);

    // assemble global stiffness matrix
    let globalk = global_k(
        points.len(),
        points[0].len(),
        material,
        &coupled_nodes,
        &mut tris,
    );

    // print the global K matrix
    println!("K =");
    print_vec2d(&globalk);
}
