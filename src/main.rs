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

// fn global_k(
//     n_nodes: usize,
//     n_freedom: usize,
//     material: (f64, f64, f64),
//     coupled_nodes: &Vec<Vec<usize>>,
//     elem_vec: &mut Vec<Triangle>,
// ) -> Vec<Vec<f64>> {
//     let mat_size: usize = n_nodes * n_freedom;
//     let mut K: Vec<Vec<f64>> = vec![vec![0.0; mat_size]; mat_size];
//
//     if coupled_nodes.len() != elem_vec.len() {
//         println!("---> Error! From assemble_global_k func.");
//         println!("     The count of elements not equal to K mat size.");
//         panic!("---> Global K failed!");
//     }
//
//     // 计算并缓存每个单元的刚度矩阵
//     let ks: Vec<_> = elem_vec.iter_mut().map(|x| x.k(material)).collect();
//
//     // 获取单个单元内的节点数目N,构造0到N的range
//     // 用于遍历单个单元的局部刚度矩阵k
//     let n_nodes: Vec<usize> = (0..coupled_nodes[0].len()).collect();
//     let loc: Vec<Vec<usize>> = full_combination(&n_nodes);
//     println!("loc = {:?}", loc);
//
//     // 整体刚度矩阵中需要修改的节点坐标对
//     // 注意！这种写法默认传进来的coupled_nodes中节点编号从1起
//     let Loc: Vec<Vec<Vec<usize>>> = coupled_nodes.iter().map(|x| full_combination(&x)).collect();
//     println!("Loc = {:?}", Loc);
//
//     for i in 0..Loc.len() {
//         for j in 0..loc.len() {
//             for k in 0..n_freedom {
//                 for l in 0..n_freedom {
//                     K[(Loc[i][j][0] - 1) * n_freedom + k][(Loc[i][j][1] - 1) * n_freedom + l] = K
//                         [(Loc[i][j][0] - 1) * n_freedom + k][(Loc[i][j][1] - 1) * n_freedom + l]
//                         + ks[i][loc[j][0] * n_freedom + k][loc[j][1] * n_freedom + l];
//                 }
//             }
//         }
//     }
//     K
// }
