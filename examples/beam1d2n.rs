use std::collections::HashMap;
use std::time::Instant;

use zhmfem::*;

fn main() {
    // set time start
    let time_start = Instant::now();

    // set material parameters
    let moment_of_inertia = 1.0 as Dtype;
    let material = (1.0 as Dtype, 0.25 as Dtype);

    let coords: Vec<Vec<Dtype>> = vec![vec![0.0, 0.0], [1.0, 0.0]];
    let cpld: Vec<Vec<usize>> = vec![vec![1, 2]];
}
