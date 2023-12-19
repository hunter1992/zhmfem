use std::time::Instant;

use zhmfem::*;

fn main() {
    let time_start = Instant::now();

    let material = (200000.0 as Dtype, 0.25 as Dtype);

    let node1 = Node2D::new(1, [-200.0, 0.0]);
    let node2 = Node2D::new(2, [0.0, 150.0]);
    let rod1: Rod2D2NNL = Rod2D2NNL::new(1, 0.1, [&node1, &node2]);

    let id = rod1.get_id();
    println!("id = {}", id);
    let l1 = rod1.length().2;
    println!("length = {}", l1);

    let disp0 = [0.0, 0.0, 0.0, 0.0];
    let b_mat = rod1.strain_matrix(&disp0);
    println!("b_mat = {:?}", b_mat);

    let km = rod1.k_m_matrix(material, &disp0);
    print_2darr("km", &km, 4.0);

    let total_time = time_start.elapsed();
    println!(">>> Total time consuming: {:?}", total_time);
}
