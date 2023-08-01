use std::collections::HashMap;
use std::time::Instant;

use zhmfem::*;

fn main() {
    // set time start
    let time_start = Instant::now();

    // set material parameters
    let moi = 0.0001186 as Dtype;
    let material = ((2.0 * 100000000000.0) as Dtype, 0.25 as Dtype);

    /*
    let coords: Vec<Vec<Dtype>> = vec![vec![0.0, 0.0], vec![1.0, 0.0]];
    let cpld: Vec<Vec<usize>> = vec![vec![1, 2]];
    let zero_disp: Vec<usize> = vec![0, 1];
    let force_index: Vec<usize> = vec![2];
    let force_value: Vec<Dtype> = vec![10.0];
    let force_data: HashMap<usize, Dtype> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();
    */

    let node1 = Node2D::new(0, [0.0, 0.0]);
    let node2 = Node2D::new(1, [5.0, 0.0]);
    let node3 = Node2D::new(1, [7.5, 0.0]);

    let beam1 = Beam1D2N::new(1, moi, 0.00665, [&node1, &node2]);
    let beam2 = Beam1D2N::new(2, moi, 0.00665, [&node2, &node3]);
    print!("{}", beam1);
    print!("{}", beam1);

    let k1 = beam1.calc_k(material);
    let k2 = beam2.calc_k(material);
    print_2darr("beam k mat", &k1);
    print_2darr("beam k mat", &k2);

    let total_time = time_start.elapsed();
    println!("\n>>> Total time consuming: {:?}", total_time);
}
