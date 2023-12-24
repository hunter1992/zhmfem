use std::time::Instant;

use zhmfem::*;

fn main() {
    let time_start = Instant::now();

    let material = (200000.0 as Dtype, 0.25 as Dtype);
    let guess_step: [Dtype; 4] = [0.0, 0.0, 0.0, 0.0];

    let node1 = Node2D::new(1, [-200.0, 0.0]);
    let node2 = Node2D::new(2, [0.0, 150.0]);
    let rod1: Rod2D2NNL = Rod2D2NNL::new(1, 100.0, [&node1, &node2]);

    let l1 = rod1.length_init();
    println!("length = {}\n", l1);

    let disp0 = [0.0, 0.0, 0.0, 0.0];
    let b_mat = rod1.strain_matrix();
    println!("b_mat = {:?}", b_mat);

    let km = rod1.k_m_stiffness_mat(material);
    print_2darr("km", &km, 4.0);

    let eg = rod1.strain_green();
    println!("strain green = {}", eg);

    let kg = rod1.k_g_stiffness_mat(material);
    print_2darr("kg", &kg, 0.0);

    let total_time = time_start.elapsed();
    println!(">>> Total time consuming: {:?}", total_time);
}
