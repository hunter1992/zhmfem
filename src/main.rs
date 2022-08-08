//extern crate nalgebra as na;

use zhmfem::*;

fn main() {
    run();
}

fn run() {
    let ee = 1.0f64;
    let nu = 0.25f64;
    let t = 1.0f64;
    let material = (ee, nu, t);

    let points: Vec<Vec<f64>> = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![1.0, 1.0],
        vec![0.0, 1.0],
    ];
    let nodes = nodes2d_vec(&points);

    let coupled_nodes: Vec<Vec<usize>> = vec![vec![1, 2, 4], vec![2, 3, 4]];
    let mut tris: Vec<Triangle> = tri2d3n_vec(&nodes, &coupled_nodes);

    let globalk = global_k(4, 2, material, &coupled_nodes, &mut tris);
    print_vec2d(&globalk);
}
