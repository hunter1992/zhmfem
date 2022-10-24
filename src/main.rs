extern crate nalgebra as na;

use na::*;
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
    for i in tris.iter_mut() {
        println!("{}", i);
        print_2darr("k", i.k(material));
    }

    // assemble global stiffness matrix
    let k_arr = global_k::<4, 2>(material, &coupled_nodes, &mut tris);

    // print the global K matrix
    print_2darr("K", &k_arr);

    // 构造nalgebra矩阵准备计算
    let k = SMatrix::<f64, 8, 8>::from(k_arr);

    // 用Gauss积分求单元的刚度矩阵

    // 构造节点位移、约束力、外力列向量
    let mut qe = vec![0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0]; // boolean
    let qe_nonzero_idx = nonzero_index(&qe);
    let fe = vector![0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    let fe_eff = fe.select_rows(qe_nonzero_idx.iter());
    let k_eff = k
        .select_columns(qe_nonzero_idx.iter())
        .select_rows(qe_nonzero_idx.iter());

    // solve the K.q = F
    let qe_unknown: Vec<f64> = k_eff.lu().solve(&fe_eff).unwrap().data.into();

    // completa qe by put qe_unknown into right position
    let _ = qe_nonzero_idx
        .iter()
        .enumerate()
        .map(|(i, &e)| qe[e] = qe_unknown[i])
        .collect::<Vec<_>>();
    print_1dvec("qe", &qe);

    tris[0].k_printer(material);
}
