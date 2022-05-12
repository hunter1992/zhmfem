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

    let elements: Vec<Vec<usize>> = vec![vec![0, 1, 3], vec![2, 3, 1]];

    let mut tri1 = Triangle::new(1, [&nodes[0], &nodes[1], &nodes[3]]);
    let mut tri2 = Triangle::new(2, [&nodes[2], &nodes[3], &nodes[1]]);

    tri1.k_printer(parameters);
    tri2.k_printer(parameters);

    //assembly(vec![tri1_k, tri2_k]);
}
