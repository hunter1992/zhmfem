/// 本例子来自曾攀《有限元分析基础教程》例题3.3.2(4)
use std::collections::HashMap;
use std::time::Instant;

use zhmfem::*;

fn main() {
    // set time start
    let time_start = Instant::now();

    // set material parameters
    let moi = 0.0001186 as Dtype;
    let cross_area = 0.00665 as Dtype;
    let material = ((2.0 * 100000000000.0) as Dtype, 0.25 as Dtype);

    let points: Vec<Vec<Dtype>> = vec![vec![0.0, 0.0], vec![5.0, 0.0], vec![7.5, 0.0]];
    let cpld = vec![vec![0, 1], vec![1, 2]];
    let zero_disp: Vec<usize> = vec![0, 1, 2];
    let force_index: Vec<usize> = vec![0, 1, 2, 3, 4, 5];
    let force_value: Vec<Dtype> = vec![-62500.0, -52083.0, -93750.0, 39062., -31250., 13021.];
    let force_data: HashMap<usize, Dtype> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    let nodes = nodes2d_vec(&points, &force_data, false);
    let mut beams: Vec<Beam1D2N> = beam1d2n_vec(moi, cross_area, &nodes, &cpld);
    print!("{}", &beams[0]);
    print!("{}", &beams[1]);

    let mut part: Part2D<Beam1D2N, 3, 2, 2> = Part2D::new(1, &nodes, &mut beams, &cpld, &material);
    part.k(material);
    part.k_printer(6.0);

    let mut eqs: LinearEqs<6> =
        LinearEqs::new(part.disps(), part.forces(), zero_disp, *part.k(material));

    eqs.lu_direct_solver();
    //eqs.gauss_seidel_iter_solver(0.000001);
    part.write_result(&eqs);

    print_1darr("qe", &part.disps(), -3.0);
    print_1darr("fe", &part.forces(), 0.0);

    let total_time = time_start.elapsed();
    println!("\n>>> Total time consuming: {:?}", total_time);
}
