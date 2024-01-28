use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::time::Instant;

use zhmfem::*;

fn main() {
    // set time start
    let time_start = Instant::now();

    let l = 0.1 as Dtype;
    let areas: [Dtype; 2] = [0.0002, 0.0001];
    let material = (2.0 as Dtype * 10000000.0 as Dtype, 0.25 as Dtype);

    const R: usize = 1;
    const C: usize = 3;
    const M: usize = 2;
    const F: usize = 1;

    let points: Vec<Vec<Dtype>> = vec![vec![0.0], vec![l], vec![2.0 * l]];
    let cpld = vec![vec![0, 1], vec![1, 2]];
    let zero_disp: Vec<usize> = vec![0];
    let force_idx: Vec<usize> = vec![2];
    let force_vlu: Vec<Dtype> = vec![10.0];
    let force_data: HashMap<usize, Dtype> =
        force_idx.into_iter().zip(force_vlu.into_iter()).collect();

    let nodes: Vec<Node1D> = nodes1d_vec(&points, &force_data);

    //let mut rod_vec: Vec<Rod1D2N> = rod1d2n_vec(section_area, &nodes, &cpld);
    let mut rod_vec: Vec<Rod1D2N> = rod1d2n_vec(&areas, &nodes, &cpld);

    let mut part1: Part1D<Rod1D2N, { R * C }, F, M> =
        Part1D::new(1, &nodes, &mut rod_vec, &cpld, &material);
    part1.k(material);
    part1.k_printer(3.0);

    let mut eqs: LinearEqs<{ R * C * F }> =
        LinearEqs::new(part1.disps(), part1.forces(), zero_disp, *part1.k(material));

    eqs.lu_direct_solver();

    part1.write_result(&eqs);

    print_1darr("qe", &part1.disps(), -3.0);
    print_1darr("fe", &part1.forces(), 0.0);

    for i in 0..rod_vec.len() {
        print!(
            "bar #{}: \nstrain:{}  stress:{}\n",
            i,
            rod_vec[i].strain([1., 2., 3.])[0],
            rod_vec[i].stress([1., 2., 3.], material)[0],
        );
    }

    /*
    let filename = "/home/zhm/Desktop/test_rod1d2n.txt";
    let file = std::fs::File::create(filename).unwrap();
    let mut writer = BufWriter::new(file);
    write!(writer, "{}", rod1.info()).expect("Write failed!");
    */

    let total_time = time_start.elapsed();
    println!(">>> Total time consuming: {:?}", total_time);
}
