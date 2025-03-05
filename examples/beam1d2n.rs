// 算例来源：曾攀《有限元初步》3.3.2节中的例4（Page59）
// 注：此例中的单位为米-牛顿

#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use std::time::Instant;

use zhmfem::*;

fn main() {
    // set timing start
    let time_start = Instant::now();

    // -------- Part 1:  Set initial parameters --------
    let cross_section_area: Vec<Dtype> = vec![0.00665, 0.00665];
    let moment_of_inertia: Vec<Dtype> = vec![0.0001186, 0.0001186];
    let material: (Dtype, Dtype) = (2.0 * 100000000000.0, 0.25); //Young's modulud & Poisson's ratio

    // Set mesh and freedom parameters
    const R: usize = 1; // rows of nodes
    const C: usize = 3; // columns of nodes
    const M: usize = 2; // num of nodes in single element
    const F: usize = 2; // num of degree freedom at single node

    //Controls the style of printing numbers in scientific notation
    const E: Dtype = 6.0;

    // Manually set coords and grouped nodes index
    let points: Vec<Vec<Dtype>> = vec![vec![0.0, 0.0], vec![5.0, 0.0], vec![7.5, 0.0]];
    let grpdnidx: Vec<Vec<usize>> = vec![vec![0, 1], vec![1, 2]];

    // Set boundary conditions and external loads
    let zero_disp_index: Vec<usize> = vec![0, 1, 2];
    let force_index: Vec<usize> = vec![0, 1, 2, 3, 4, 5];
    let force_value: Vec<Dtype> = vec![-62500.0, -52083.0, -93750.0, 39062.0, -31250.0, 13021.0];
    let force_data: HashMap<usize, Dtype> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // -------- Part 2:  Construct nodes, elements and parts --------
    // Construct 2D nodes vector
    let nodes = nodes2d_vec(&points, &force_data);

    // Construct Beam1D2N elements vector
    let mut rods = beam1d2n_vec(
        &nodes,
        &grpdnidx,
        &moment_of_inertia,
        &cross_section_area,
        &material,
    );
    let element_type: &str = "Beam1D2N_";

    // Construct 2D part & assembly global stiffness matrix
    let mut part: Part2D<'_, Beam1D2N<'_>, { R * C }, F, M> =
        Part2D::new(1, &nodes, &mut rods, &grpdnidx);
    let parallel_or_singllel = "singllel";
    part.k_printer(parallel_or_singllel, E);

    // -------- Part 3:  Solve the problem --------
    // construct solver and solve the case
    let mut eqs: LinearEqs<{ R * C * F }> = LinearEqs::new(
        part.nodes_displacement(),
        part.nodes_force(),
        zero_disp_index,
        *part.k(parallel_or_singllel),
    );

    // 1) solve the linear equations of static system using direct method.
    eqs.lu_direct_solver(); //LU decomposition method
    let output_file = "LU.txt";

    // 2) solve the linear equations of static system using iter method.
    //eqs.gauss_seidel_iter_solver(0.001);
    //let output_file = "G-S.txt";

    let calc_time: std::time::Duration = eqs.solver_time_consuming.unwrap();

    // write the displacement and force result into Node2D's field
    part.write_result(&eqs);

    // -------- Part 4:  Print all kinds of result --------
    print_1darr("qe", &part.nodes_displacement(), 0.0, "v");
    print_1darr("fe", &part.nodes_force(), 4.0, "v");

    println!("\n>>> System energy:");
    let strain_energy: Dtype =
        strain_energy(*part.k(parallel_or_singllel), part.nodes_displacement());
    let external_force_work: Dtype =
        external_force_work(part.nodes_force(), part.nodes_displacement());
    let potential_energy: Dtype = potential_energy(
        *part.k(parallel_or_singllel),
        part.nodes_force(),
        part.nodes_displacement(),
    );
    println!("\tE_d: {:-9.6} (deform energy)", strain_energy);
    println!("\tW_f: {:-9.6} (exforce works)", external_force_work);
    println!("\tE_p: {:-9.6} (potential energy)", potential_energy);

    for elem in part.elems.iter() {
        elem.k_printer(E);
        elem.print_strain([2.5, 0.317, 0.0]);
        elem.print_stress([2.5, 0.317, 0.0]);
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
