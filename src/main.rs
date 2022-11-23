extern crate nalgebra as na;

use na::*;
use std::collections::HashMap;
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

    // input the node coordinates
    let coords: Vec<Vec<f64>> = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![1.0, 1.0],
        vec![0.0, 1.0],
    ];
    let disp_0_idx: Vec<usize> = vec![0, 1, 6];
    let force_index: Vec<usize> = vec![2, 4];
    let force_value: Vec<f64> = vec![-1.0, 1.0];
    let force_data: HashMap<usize, f64> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // transform points into nodes
    let nodes = nodes2d_vec(&coords, &disp_0_idx, force_data);

    // list nodes ids in one element
    let cpld: Vec<Vec<usize>> = vec![vec![1, 2, 4], vec![3, 4, 2]];

    // construct element by coupled nodes
    let mut tris: Vec<Tri2D3N> = tri2d3n_vec(&nodes, &cpld);
    for i in tris.iter_mut() {
        println!("{}", i);
        i.k_printer(material);
    }

    // assemble global stiffness matrix
    let mut p1: Part2D<Tri2D3N, 4, 2, 3> = Part2D::new(1, tris, cpld, material);

    // print the global K matrix
    print_2darr("K", p1.k());

    // 构造nalgebra矩阵准备计算
    let k = SMatrix::<f64, 8, 8>::from(*p1.k());

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

    // complete qe by put qe_unknown into right position
    let _ = qe_nonzero_idx
        .iter()
        .enumerate()
        .map(|(i, &e)| qe[e] = qe_unknown[i])
        .collect::<Vec<_>>();
    print_1dvec("qe", &qe);
}
