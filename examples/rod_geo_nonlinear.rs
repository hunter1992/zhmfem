use std::collections::HashMap;
use std::time::Instant;

use zhmfem::*;

fn main() {
    const N: usize = 3;
    const M: usize = 2;
    const F: usize = 2;

    let time_start = Instant::now();

    let section_area: [Dtype; 2] = [100.0; 2];
    let material = (200000.0 as Dtype, 0.25 as Dtype);

    let points: Vec<Vec<Dtype>> = vec![vec![-200.0, 0.0], vec![0.0, 150.0], vec![200.0, 0.0]];
    let cpld = vec![vec![0, 1], vec![2, 1]];

    let zero_disp: Vec<usize> = vec![0, 1, 2, 3, 4, 5];
    let force_idx: Vec<usize> = vec![2, 3];
    let force_vlu: Vec<Dtype> = vec![0.0, 1000000.0];
    let force_data: HashMap<usize, Dtype> =
        force_idx.into_iter().zip(force_vlu.into_iter()).collect();

    let nodes: Vec<Node2D> = nodes2d_vec(&points, &force_data, true);

    let mut rods: Vec<Rod2D2NNL> = rod2d2n_nonlinear_vec(&section_area, &nodes, &cpld);
    for rod in rods.iter_mut() {
        rod.write_node_inner_force(material);
    }

    let mut part1: Part2D<Rod2D2NNL, N, F, M> = Part2D::new(1, &nodes, &mut rods, &cpld, &material);
    let global_kt = part1.k(material);
    print_2darr("Kt", global_kt, 4.0);
    print_1darr("inner force", &part1.forces(), 0.0);
    let ext_force = part1.nodes_external_forces(&force_data);
    print_1darr("External force", &ext_force, 0.0);

    let total_time = time_start.elapsed();
    println!(">>> Total time consuming: {:?}", total_time);
}
