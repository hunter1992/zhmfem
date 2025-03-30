#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

/** This example comes from:
 *  https://cloud.tencent.com/developer/article/1086452
 */
use std::time::Instant;

use zhmfem::*;

fn main() {
    // set timing start
    let time_start = Instant::now();

    // -------- Part 1:  Set initial parameters --------
    let thick: Dtype = 10.0; //Thickness of the plate
    let material: (Dtype, Dtype) = (206000.0, 0.3); //Young's modulud & Poisson's ratio

    // Set mesh and freedom parameters
    const R: usize = 2; // rows of nodes
    const C: usize = 2; // columns of nodes
    const M: usize = 3; // num of nodes in single element
    const F: usize = 2; // num of degree freedom at single node
    const CPU_CORES: usize = 2;

    //Controls the style of printing numbers in scientific notation
    const E: Dtype = 6.0;

    // Manually set coords and grouped nodes index
    const L: Dtype = 300.;
    const H: Dtype = 200.;
    let points: Vec<Vec<Dtype>> = vec![vec![L, 0.0], vec![L, H], vec![0.0, H], vec![0.0, 0.0]];
    let grpdnidx: Vec<Vec<usize>> = vec![vec![0, 1, 3], vec![2, 3, 1]];

    // Set boundary conditions and external loads
    let zero_disp_index: Vec<usize> = vec![1, 4, 5, 6, 7];
    let force_index: Vec<usize> = vec![3];
    let force_value: Vec<Dtype> = vec![-500000.0];
    let force_data: HashMap<usize, Dtype> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // -------- Part 2:  Construct nodes, elements and parts --------
    // Construct 2D nodes vector
    let nodes = nodes2d_vec(&points, &force_data);

    // Construct Tri2D3N elements vector
    let mut triangles = tri2d3n_vec(thick, &nodes, &grpdnidx, &material);
    //let element_type: &str = "Tri2D3N_test_";

    // Construct 2D part & assembly global stiffness matrix
    let mut part: Part2D<'_, Tri2D3N<'_>, { R * C }, F, M> =
        Part2D::new(1, &nodes, &mut triangles, &grpdnidx);
    let parallel_or_singllel: &str = "s";
    part.k_printer(parallel_or_singllel, CPU_CORES, E);

    // -------- Part 3:  Solve the problem --------
    // construct solver and solve the case
    let mut eqs: LinearEqs<{ R * C * F }> = LinearEqs::new(
        part.nodes_displacement(),
        part.nodes_force(),
        zero_disp_index,
        *part.k(parallel_or_singllel, CPU_CORES),
    );

    // 1) solve the linear equations of static system using direct method.
    eqs.lu_direct_solver(); //LU decomposition method
                            //let output_file = "LU.txt";

    // 2) solve the linear equations of static system using iter method.
    //eqs.gauss_seidel_iter_solver(0.001);
    //let output_file = "G-S.txt";

    // 3) or you can solve the problem with a more concise call:
    // eqs.solve("lu", 0.001); // or eqs.solve("gs", 0.001);

    //let calc_time: std::time::Duration = eqs.solver_time_consuming.unwrap();

    // write the displacement and force result into Node2D's field
    part.write_result(&eqs);

    // -------- Part 4:  Print all kinds of result --------
    print_1darr("qe", &part.nodes_displacement(), 0., "v");
    print_1darr("fe", &part.nodes_force(), E, "v");

    println!("\n>>> System energy:");
    let strain_energy: Dtype = strain_energy(
        *part.k(parallel_or_singllel, CPU_CORES),
        part.nodes_displacement(),
    );
    let external_force_work: Dtype =
        external_force_work(part.nodes_force(), part.nodes_displacement());
    let potential_energy: Dtype = potential_energy(
        *part.k(parallel_or_singllel, CPU_CORES),
        part.nodes_force(),
        part.nodes_displacement(),
    );
    println!("\tE_d: {:-9.6} (deform energy)", strain_energy);
    println!("\tW_f: {:-9.6} (exforce works)", external_force_work);
    println!("\tE_p: {:-9.6} (potential energy)", potential_energy);

    for elem in part.elems.iter() {
        elem.k_printer(E);
        elem.print_strain();
        elem.print_stress();
    }

    // -------- Part 5:  Write clac result into txt file --------
    /*
    let output_path = "/home/zhm/Documents/Scripts/Rust/zhmfem/results/";
    let output = format!("{output_path}{element_type}{output_file}");
    part.txt_writer(
        &output,
        calc_time,
        E,
        (strain_energy, external_force_work, potential_energy),
    )
    .expect(">>> !!! Failed to output text result file !!!");
    */

    let total_time = time_start.elapsed();
    println!("\n>>> Total time consuming: {:?}", total_time);
}
