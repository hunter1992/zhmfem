/// 计算悬臂梁在静态载荷作用下的小变形挠度
/// 这个例子来源于：
/// 张大羽, 罗建君, 等. 基于旋转场曲率的二维剪切梁单元建模[J]. 物理学报, 66, 114501(2017).
/// DOI: 10.7498/aps.66.114501
use std::collections::HashMap;
use std::time::Instant;

use zhmfem::*;

fn main() {
    // set time start
    let time_start = Instant::now();

    // set material parameters
    let moi = 0.00104167 as Dtype;
    let cross_area = 0.05 as Dtype;
    let material = ((2.07 * 100000000000.0) as Dtype, 0.3 as Dtype);

    // freedom parameters
    const R: usize = 1;
    const C: usize = 2;
    const M: usize = 2;
    const F: usize = 2;

    // coords, cpld nodes and load force data
    let points: Vec<Vec<Dtype>> = vec![vec![0.0, 0.0], vec![2.0, 0.0]];
    let cpld = vec![vec![0, 1]];
    let zero_disp: Vec<usize> = vec![0, 1];
    let force_index: Vec<usize> = vec![2];
    let force_value: Vec<Dtype> = vec![-62500.0];
    let force_data: HashMap<usize, Dtype> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // construct nodes and elements
    let nodes = nodes2d_vec(&points, &force_data, false);
    let mut beam_vec: Vec<Beam1D2N> = beam1d2n_vec(moi, cross_area, &nodes, &cpld);
    print!("{}", &beam_vec[0]);

    // construct part
    let mut beam_part: Part2D<Beam1D2N, { R * C }, F, M> =
        Part2D::new(1, &nodes, &mut beam_vec, &cpld, &material);
    beam_part.k(material);
    beam_part.k_printer(0.0);

    // construct solver
    let mut eqs: LinearEqs<{ R * C * F }> = LinearEqs::new(
        beam_part.disps(),
        beam_part.forces(),
        zero_disp,
        *beam_part.k(material),
    );

    // solving the problem and write result into elements field
    //eqs.lu_direct_solver();
    eqs.gauss_seidel_iter_solver(0.0001);
    beam_part.write_result(&eqs);

    print_1darr("qe", &beam_part.disps(), -4.0);
    print_1darr("fe", &beam_part.forces(), 0.0);

    let total_time = time_start.elapsed();
    println!("\n>>> Total time consuming: {:?}", total_time);
}
