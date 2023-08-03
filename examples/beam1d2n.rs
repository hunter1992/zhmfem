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
    let cpld = vec![vec![1, 2], vec![2, 3]];
    let zero_disp: Vec<usize> = vec![0, 1, 2];
    let force_index: Vec<usize> = vec![3, 4, 5];
    let force_value: Vec<Dtype> = vec![39062., -31250., 13021.];
    let force_data: HashMap<usize, Dtype> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    let nodes = nodes2d_vec(&points, &zero_disp, &force_data);
    let mut beam_vec: Vec<Beam1D2N> = beam1d2n_vec(moi, cross_area, &nodes, &cpld);
    print!("{}", &beam_vec[0]);
    print!("{}", &beam_vec[1]);

    let mut beam_part: Part2D<Beam1D2N, 3, 2, 2> = Part2D::new(1, &nodes, &mut beam_vec, &cpld);
    beam_part.k(material);
    beam_part.k_printer(6.0);

    let mut eqs: LinearEqs<6> = LinearEqs::new(
        beam_part.disps(),
        beam_part.forces(),
        *beam_part.k(material),
    );

    eqs.lu_direct_solver();
    beam_part.write_result(&eqs);

    print_1darr("qe", &beam_part.disps());
    print_1darr("fe", &beam_part.forces());

    let total_time = time_start.elapsed();
    println!("\n>>> Total time consuming: {:?}", total_time);
}
