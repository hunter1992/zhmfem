#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use std::time::Instant;

use zhmfem::*;

fn main() {
    // set timing start
    let time_start = Instant::now();

    // -------- Part 1:  Set initial parameters --------
    let cross_section_area: Vec<Dtype> = vec![100.0, 100.0, 100.0, 100.0];
    let material: (Dtype, Dtype) = (295000.0, 0.25); //Young's modulud & Poisson's ratio

    // Set mesh and freedom parameters
    const R: usize = 2; // rows of nodes
    const C: usize = 2; // columns of nodes
    const M: usize = 2; // num of nodes in single element
    const F: usize = 2; // num of degree freedom at single node

    //Controls the style of printing numbers in scientific notation
    const E: Dtype = 4.0;

    // Manually set coords and grouped nodes index
    let points: Vec<Vec<Dtype>> = vec![
        vec![0.0, 0.0],
        vec![400.0, 0.0],
        vec![400.0, 300.0],
        vec![0.0, 300.0],
    ];
    let grpdnidx: Vec<Vec<usize>> = vec![vec![0, 1], vec![2, 1], vec![0, 2], vec![3, 2]];

    // Set boundary conditions and external loads
    let zero_disp_index: Vec<usize> = vec![0, 1, 3, 6, 7];
    let force_index: Vec<usize> = vec![2, 5];
    let force_value: Vec<Dtype> = vec![20000.0, -25000.0];
    let force_data: HashMap<usize, Dtype> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // -------- Part 2:  Construct nodes, elements and parts --------
    // Construct 2D nodes vector
    let nodes = nodes2d_vec(&points, &force_data);

    // Construct Rod2D2N elements vector
    let mut rods = rod2d2n_vec(&nodes, &grpdnidx, &cross_section_area, &material);
    let element_type: &str = "Rod2D2N_";

    // Construct 2D part & assembly global stiffness matrix
    let mut part: Part2D<'_, Rod2D2N<'_>, { R * C }, F, M> =
        Part2D::new(1, &nodes, &mut rods, &grpdnidx);
    part.k_printer(E);

    // -------- Part 3:  Solve the problem --------
    // construct solver and solve the case
    let mut eqs: LinearEqs<{ R * C * F }> = LinearEqs::new(
        part.nodes_displacement(),
        part.nodes_force(),
        zero_disp_index,
        *part.k(),
    );

    // 1) solve the linear equations of static system using direct method.
    eqs.lu_direct_solver(); //LU decomposition method
    let output_file = "LU.txt";

    // 2) solve the linear equations of static system using iter method.
    //eqs.gauss_seidel_iter_solver(0.001);
    //let output_file = "G-S.txt";

    let calc_time: std::time::Duration = eqs.solver_calc_time.unwrap();

    // write the displacement and force result into Node2D's field
    part.write_result(&eqs);

    // -------- Part 4:  Print all kinds of result --------
    print_1darr("qe", &part.nodes_displacement(), 0.0);
    print_1darr("fe", &part.nodes_force(), E);

    println!("\n>>> System energy:");
    let strain_energy: Dtype = strain_energy(*part.k(), part.nodes_displacement());
    let external_force_work: Dtype =
        external_force_work(part.nodes_force(), part.nodes_displacement());
    let potential_energy: Dtype =
        potential_energy(*part.k(), part.nodes_force(), part.nodes_displacement());
    println!("\tE_d: {:-9.6} (deform energy)", strain_energy);
    println!("\tW_f: {:-9.6} (exforce works)", external_force_work);
    println!("\tE_p: {:-9.6} (potential energy)", potential_energy);

    for elem in part.elems.iter() {
        elem.k_printer(E);
        elem.print_strain();
        elem.print_stress();
    }

    // -------- Part 5:  Write clac result into txt file --------
    let output_path = "/home/zhm/Documents/Scripts/Rust/zhmfem/results/";
    let output = format!("{output_path}{element_type}{output_file}");
    part.txt_writer(
        &output,
        calc_time,
        E,
        (strain_energy, external_force_work, potential_energy),
    )
    .expect(">>> !!! Failed to output text result file !!!");

    let total_time = time_start.elapsed();
    println!("\n>>> Total time consuming: {:?}", total_time);
}
