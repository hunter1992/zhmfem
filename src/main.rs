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

    let K = global_k(4, 2, material, &coupled_nodes, &mut tris);
    //vec2d_printer(&kk);

    //let mut tri1 = Triangle::new(1, [&nodes[0], &nodes[1], &nodes[3]]);
    //let mut tri2 = Triangle::new(2, [&nodes[2], &nodes[3], &nodes[1]]);
}

fn global_k(
    n_nodes: usize,
    n_freedom: usize,
    material: (f64, f64, f64),
    coupled_nodes: &Vec<Vec<usize>>,
    elem_vec: &mut Vec<Triangle>,
) -> Vec<Vec<f64>> {
    let mat_size: usize = n_nodes * n_freedom;
    let mut K: Vec<Vec<f64>> = vec![vec![0.0; mat_size]; mat_size];

    if coupled_nodes.len() != elem_vec.len() {
        println!("---> Error! From assemble_global_k func.");
        println!("     The count of elements not equal to K mat size.");
        panic!("---> Global K failed!");
    }

    // 计算并缓存每个单元的刚度矩阵
    elem_vec.iter_mut().map(|x| x.k(material));

    // 获取单个单元内的节点数目N,构造0到N的range
    // 用于遍历单个单元的局部刚度矩阵k
    let n_nodes: Vec<usize> = (0..coupled_nodes[0].len()).collect();
    let loc: Vec<Vec<usize>> = full_combination(&n_nodes);

    // 整体刚度矩阵中需要修改的节点坐标对
    // 注意！这种写法默认传进来的coupled_nodes中节点编号从1起
    let Loc: Vec<Vec<Vec<usize>>> = coupled_nodes.iter().map(|x| full_combination(&x)).collect();

    for i in 0..Loc.len(){
        for j in 0..loc.len(){
           K[] 
        }
    }
    K
}

//fn assembly(aim_mat: &mut Vec<Vec<f64>>, aim_num: &f64, loc: &Vec<(usize, usize)>) {
//    aim
//}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
