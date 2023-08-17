use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::time::Instant;

use zhmfem::*;

fn main() {
    let time_start = Instant::now();

    let l = 0.1 as Dtype;
    let area1 = 0.02 as Dtype;
    let area2 = 0.03 as Dtype;
    let area3 = 0.06 as Dtype;
    let material = (2.0 as Dtype * 100000.0 as Dtype, 0.25 as Dtype);

    const R: usize = 1;
    const C: usize = 4;
    const M: usize = 2;
    const F: usize = 1;

    let points: Vec<Vec<Dtype>> = vec![vec![0.0], vec![l], vec![2.0 * l], vec![3.0 * l]];
    let cpld = vec![vec![0, 1], vec![1, 2], vec![2, 3]];
    let zero_disp: Vec<usize> = vec![3];
    let force_idx: Vec<usize> = vec![0, 2];
    let force_vlu: Vec<Dtype> = vec![-100.0, 50.0];
    let force_data: HashMap<usize, Dtype> =
        force_idx.into_iter().zip(force_vlu.into_iter()).collect();

    let nodes: Vec<Node1D> = nodes1d_vec(&points, &zero_disp, &force_data);

    //let mut rod_vec: Vec<Rod1D2N> = rod1d2n_vec(section_area, &nodes, &cpld);
    let mut rod_vec: Vec<Rod1D2N> = vec![
        Rod1D2N::new(1, area1, [&nodes[0], &nodes[1]]),
        Rod1D2N::new(2, area2, [&nodes[1], &nodes[2]]),
        Rod1D2N::new(3, area3, [&nodes[2], &nodes[3]]),
    ];

    let mut part1: Part1D<Rod1D2N, { R * C }, F, M> = Part1D::new(1, &nodes, &mut rod_vec, &cpld);
    part1.k(material);
    part1.k_printer(0.0);

    let mut eqs: LinearEqs<{ R * C * F }> =
        LinearEqs::new(part1.disps(), part1.forces(), *part1.k(material));

    eqs.lu_direct_solver();

    part1.write_result(&eqs);

    print_1darr("qe", &part1.disps(), -3.0);
    print_1darr("fe", &part1.forces(), 0.0);

    /*
    let filename = "/home/zhm/Desktop/test_rod1d2n.txt";
    let file = std::fs::File::create(filename).unwrap();
    let mut writer = BufWriter::new(file);
    write!(writer, "{}", rod1.info()).expect("Write failed!");
    */

    let total_time = time_start.elapsed();
    println!(">>> Total time consuming: {:?}", total_time);
}
