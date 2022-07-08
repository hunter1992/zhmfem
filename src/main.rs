//extern crate nalgebra as na;

//use na::Vector3;
use zhmfem::*;

fn main() {
    run();
}

fn run() {
    let ee = 1.0f64;
    let nu = 0.25f64;
    let t = 1.0f64;
    let parameters = (ee, nu, t);

    let points: Vec<Vec<f64>> = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![1.0, 1.0],
        vec![0.0, 1.0],
    ];
    let nodes = nodes2d_vec(&points);

    let mut tri1 = Triangle::new(1, [&nodes[0], &nodes[1], &nodes[3]]);
    let mut tri2 = Triangle::new(2, [&nodes[2], &nodes[3], &nodes[1]]);

    k1 = tri1.k(parameters);
    k2 = tri2.k(parameters);

    let K = assembly(k1, k2);
}

fn assembly_tri2D3N(
    n_nodes: usize,
    n_freedom: usize,
    coupled_nodes: Vec<Vec<usize>>,
    ks: Vec<Vec<&[[f64; 6]; 6]>>,
) {
    let k
}
