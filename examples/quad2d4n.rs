#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use std::time::Instant;

use zhmfem::*;

fn main() {
    // Set time start
    let time_start = Instant::now();

    // ------ Part 1: Set initial parameters ------
    const E: Dtype = 0.0; // Exponent in scientific notation to base 10
    const CPU_CORES: usize = 2;

    let calc_method: &str = "lu"; // "lu" for LU decomposition algorithm or "gs" for gauss-seidel iteration method
    let calc_accuracy: Dtype = 0.001; // Calculation accuracy of iterative algorithm

    let parallel_or_singllel: &str = "s"; // "s" or "p"

    let thick: Dtype = 0.1; //Thickness of the plate
    let force: Dtype = 10000.0;
    let material: (Dtype, Dtype) = (30.0, 0.3); //Young's modulud & Poisson's ratio

    // Number of nodes and freedom
    const R: usize = 2; // rows of nodes
    const C: usize = 2; // columns of nodes
    const M: usize = 4; // node num in single element
    const F: usize = 2; // freedom num in single node

    // Manually set coords and grouped nodes index
    let points: Vec<Vec<Dtype>> = vec![
        vec![1.0, 0.0],
        vec![2.0, 0.0],
        vec![2.25, 1.5],
        vec![1.25, 1.0],
    ];
    let grpdnidx: Vec<Vec<usize>> = vec![vec![0, 1, 2, 3]];

    // Set boundary conditions and mesh it
    let zero_disp_index: Vec<usize> = vec![0, 1, 2, 3];
    let force_index: Vec<usize> = vec![4];
    let force_value: Vec<Dtype> = vec![force];
    let force_data: HashMap<usize, Dtype> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // -------- Part 2:  Construct nodes, elements and parts --------
    // construct 2D nodes vector
    let nodes = nodes2d_vec(&points, &force_data);

    // Construct Quad2D4N elements vector
    let mut quads = quad2d4n_vec(thick, &nodes, &grpdnidx, &material);
    let element_type: &str = "Quad2D4N_";

    // Construct 2D part & assembly global stiffness matrix
    let mut part: Part2D<'_, Quad2D4N<'_>, { R * C }, F, M> =
        Part2D::new(1, &nodes, &mut quads, &grpdnidx);
    let parallel_or_singllel = "singllel";
    part.k_printer(parallel_or_singllel, CPU_CORES, E);

    // -------- Part 3:  Solve the problem --------
    // construct solver and solve the case
    let mut eqs: LinearEqs<{ R * C * F }> = LinearEqs::new(
        part.nodes_displacement(),
        part.nodes_force(),
        zero_disp_index,
        part.k(parallel_or_singllel, CPU_CORES).clone(),
    );

    // solve method:
    //     "lu" --- LU decomposition method
    //     "gs" --- gaussr seidel iter method
    //
    //eqs.solve("gs", 0.001);
    //let output_file = "G-S.txt";
    eqs.solve("lu", 0.001);
    let output_file = "LU.txt";

    // write the displacement and force result into nodes' field
    part.write_result(&eqs);

    let calc_time: std::time::Duration = eqs.solver_time_consuming.unwrap();

    // -------- Part 4:  Print all kinds of result --------
    print_1darr("qe", &part.nodes_displacement(), 0.0, "v");
    print_1darr("fe", &part.nodes_force(), E, "v");

    println!("\n>>> System energy:");
    let strain_energy: Dtype = strain_energy(
        part.k(parallel_or_singllel, CPU_CORES).clone(),
        part.nodes_displacement(),
    );
    let external_force_work: Dtype =
        external_force_work(part.nodes_force(), part.nodes_displacement());
    let potential_energy: Dtype = potential_energy(
        part.k(parallel_or_singllel, CPU_CORES).clone(),
        part.nodes_force(),
        part.nodes_displacement(),
    );
    println!("\tE_d: {:-9.6} (deform energy)", strain_energy);
    println!("\tW_f: {:-9.6} (exforce works)", external_force_work);
    println!("\tE_p: {:-9.6} (potential energy)", potential_energy);

    for elem in part.elems.iter() {
        elem.k_printer(E);
        elem.print_strain([0., 0., 0.]);
        elem.print_stress([0., 0., 0.]);
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
