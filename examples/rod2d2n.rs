use std::collections::HashMap;
//use std::io::{BufWriter, Write};
use std::time::Instant;

use zhmfem::*;

fn main() {
    let time_start = Instant::now();

    let section_area: [Dtype; 4] = [100.0; 4];
    let material = (295000.0 as Dtype, 0.25 as Dtype);

    const R: usize = 2;
    const C: usize = 2;
    const M: usize = 2;
    const F: usize = 2;

    let points: Vec<Vec<Dtype>> = vec![
        vec![0.0, 0.0],
        vec![400.0, 0.0],
        vec![400.0, 300.0],
        vec![0.0, 300.0],
    ];
    let cpld = vec![vec![0, 1], vec![2, 1], vec![0, 2], vec![3, 2]];
    let zero_disp: Vec<usize> = vec![0, 1, 3, 6, 7];
    let force_idx: Vec<usize> = vec![2, 5];
    let force_vlu: Vec<Dtype> = vec![20000.0, -25000.0];
    let force_data: HashMap<usize, Dtype> =
        force_idx.into_iter().zip(force_vlu.into_iter()).collect();

    let nodes: Vec<Node2D> = nodes2d_vec(&points, &zero_disp, &force_data);

    let mut rod_vec: Vec<Rod2D2N> = rod2d2n_vec(&section_area, &nodes, &cpld);

    let mut part1: Part2D<Rod2D2N, { R * C }, F, M> = Part2D::new(1, &nodes, &mut rod_vec, &cpld);
    part1.k(material);
    part1.k_printer(3.0);

    let mut eqs: LinearEqs<{ R * C * F }> =
        LinearEqs::new(part1.disps(), part1.forces(), *part1.k(material));

    eqs.lu_direct_solver();

    part1.write_result(&eqs);

    print_1darr("qe", &part1.disps(), 0.0);
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
