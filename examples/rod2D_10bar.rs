/// 算例来源：左文杰《一维有限单元法》7.4节中的十杆桁架结构
/// 注：此例中的单位为英寸-磅
use std::collections::HashMap;
//use std::io::{BufWriter, Write};
use std::time::Instant;

use zhmfem::*;

fn main() {
    // set time start
    let time_start = Instant::now();

    let section_area: [Dtype; 10] = [10.0; 10];
    let material = (10000000.0 as Dtype, 0.25 as Dtype);

    const R: usize = 2;
    const C: usize = 3;
    const M: usize = 2;
    const F: usize = 2;

    let l: Dtype = 360.0;
    let points: Vec<Vec<Dtype>> = vec![
        vec![0.0, l],
        vec![l, l],
        vec![2.0 * l, l],
        vec![2.0 * l, 0.0],
        vec![l, 0.0],
        vec![0.0, 0.0],
    ];
    let cpld = vec![
        vec![0, 1],
        vec![5, 4],
        vec![1, 2],
        vec![4, 3],
        vec![0, 4],
        vec![5, 1],
        vec![1, 3],
        vec![4, 2],
        vec![1, 4],
        vec![2, 3],
    ];
    let zero_disp: Vec<usize> = vec![0, 1, 10, 11];
    let force_idx: Vec<usize> = vec![2, 3, 4, 5, 6, 7, 8, 9];
    let force_vlu: Vec<Dtype> = vec![0., 0., 0., 0., 0., -100000.0, 0., -100000.0];
    let force_data: HashMap<usize, Dtype> =
        force_idx.into_iter().zip(force_vlu.into_iter()).collect();

    let nodes: Vec<Node2D> = nodes2d_vec(&points, &zero_disp, &force_data);

    let mut rod_vec: Vec<Rod2D2N> = rod2d2n_vec(&section_area, &nodes, &cpld);

    let mut part1: Part2D<Rod2D2N, { R * C }, F, M> = Part2D::new(1, &nodes, &mut rod_vec, &cpld);
    part1.k(material);
    part1.k_printer(4.0);

    let mut eqs: LinearEqs<{ R * C * F }> =
        LinearEqs::new(part1.disps(), part1.forces(), *part1.k(material));

    eqs.lu_direct_solver();
    //eqs.gauss_seidel_iter_solver(0.00001); // 耗时:1.889097ms,迭代:143次,误差:0.000009

    part1.write_result(&eqs);

    print_1darr("qe", &part1.disps(), 0.0);
    print_1darr("fe", &part1.forces(), 4.0);

    /*
    for rod in rod_vec.iter() {
        print_1dvec("stress", &rod.stress([1., 1., 1.], material), 4.);
    }


    let filename = "/home/zhm/Desktop/test_rod1d2n.txt";
    let file = std::fs::File::create(filename).unwrap();
    let mut writer = BufWriter::new(file);
    write!(writer, "{}", rod1.info()).expect("Write failed!");
    */

    let total_time = time_start.elapsed();
    println!(">>> Total time consuming: {:?}", total_time);
}
