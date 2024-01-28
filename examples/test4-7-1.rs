use std::collections::HashMap;
use std::time::Instant;

use zhmfem::*;

fn main() {
    // set timing start
    let time_start = Instant::now();

    // set material parameters
    let thick: Dtype = 0.1;
    let ex_force: Dtype = 100000.0;
    let material: (Dtype, Dtype) = (10000000.0, 0.33); //Young's modulud & Poisson's ratio

    const R: usize = 2;
    const C: usize = 2;
    const F: usize = 2;
    const M: usize = 3;

    let points: Vec<Vec<Dtype>> = vec![
        vec![2.0, 1.0],
        vec![2.0, 0.0],
        vec![0.0, 1.0],
        vec![0.0, 0.0],
    ];
    let cpld: Vec<Vec<usize>> = vec![vec![3, 1, 2], vec![1, 0, 2]];

    // set boundary conditions and loads
    let zero_disp: Vec<usize> = vec![4, 5, 6, 7];
    let force_index: Vec<usize> = vec![1, 3];
    let force_value: Vec<Dtype> = vec![-0.5 * ex_force, -0.5 * ex_force];
    let force_data: HashMap<usize, Dtype> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // transform points into nodes
    let nodes: Vec<Node2D> = nodes2d_vec(&points, &force_data, false);

    // construct elements by coupled nodes
    let mut tri_vec: Vec<Tri2D3N> = tri2d3n_vec(thick, &nodes, &cpld);

    // assemble global stiffness matrix
    let mut part1: Part2D<Tri2D3N, { R * C }, F, M> =
        Part2D::new(1, &nodes, &mut tri_vec, &cpld, &material);
    //println!("");
    part1.k(material);
    part1.k_printer(6.0);

    // construct solver and solve the case
    let mut eqs: LinearEqs<{ R * C * F }> =
        LinearEqs::new(part1.disps(), part1.forces(), zero_disp, *part1.k(material));

    // 1) solve the linear equations of static system using direct method.
    eqs.lu_direct_solver();

    // 2) solve the linear equations of static system using iter method.
    //eqs.gauss_seidel_iter_solver(0.001);

    part1.write_result(&eqs);

    print_1darr("qe", &part1.disps(), 0.0);
    print_1darr("fe", &part1.forces(), 0.0);

    for tri in part1.elems.iter() {
        tri.print_strain();
        tri.print_stress(material);
    }

    println!("\n>>> System energy:");
    println!("\tE_d: {:-9.6} (deform energy)", part1.strain_energy());
    println!("\tW_f: {:-9.6} (exforce works)", part1.force_work());
    println!(
        "\tE_p: {:-9.6} (potential energy)",
        part1.potential_energy()
    );

    // write clac result into txt file
    //part1.write_txt_file(material, "/home/zhm/Desktop/tri.txt");

    let total_time = time_start.elapsed();
    println!("\n>>> Total time consuming: {:?}", total_time);
}
