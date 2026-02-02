// 这个例子来源于曾攀《有限元基础教程》的例题3.2.5(1)

use std::time::Instant;

use zhmfem::*;

fn main() {
    // set timing start
    let time_start = Instant::now();

    // "lu"       for LU       decomposition algorithm or
    // "cholesky" for Cholesky decomposition algorithm or
    // "pardiso"  for calling  Panua Tech's PARDISO
    // "gs"       for gauss-seidel iteration algorithm
    let calc_method: &str = "cholesky";
    let calc_accuracy: Dtype = 0.0001;
    let parallel_or_singllel: &str = "s";

    // -------- Part 1:  Set initial parameters --------
    let cross_section_area: Vec<Dtype> = vec![100.0, 100.0, 100.0, 100.0];
    let material: [Dtype; 2] = [295000.0, 0.25]; //Young's modulud & Poisson's ratio

    // Set mesh and freedom parameters
    const R: usize = 2; // rows of nodes
    const C: usize = 2; // columns of nodes
    const M: usize = 2; // num of nodes in single element
    const F: usize = 2; // num of degree freedom at single node
    const CPU_CORES: usize = 2;

    //Controls the style of printing numbers in scientific notation
    const E: Dtype = 4.0;

    // Manually set coords and grouped nodes index
    let points: Vec<[Dtype; 2]> = vec![[0.0, 0.0], [400.0, 0.0], [400.0, 300.0], [0.0, 300.0]];
    let grpdnidx: Vec<Vec<usize>> = vec![vec![0, 1], vec![2, 1], vec![0, 2], vec![3, 2]];

    // Set boundary conditions and external loads
    let zero_disp_index: Vec<usize> = vec![0, 1, 3, 6, 7];
    let force_index: Vec<usize> = vec![2, 5];
    let force_value: Vec<Dtype> = vec![20000.0, -25000.0];

    // -------- Part 2:  Construct nodes, elements and parts --------
    // Construct 2D nodes vector
    let nodes: Vec<Node2D> = nodes2d_vec(&points, &force_index, &force_value);

    // Construct Rod2D2N elements vector
    let mut rods = rod2d2n_vec(&nodes, &grpdnidx, &cross_section_area, material);

    // Construct 2D part & assembly global stiffness matrix
    let mut part: Part2D<'_, Rod2D2N<'_>, { R * C }, F, M> =
        Part2D::new(1, &nodes, &mut rods, &grpdnidx);
    part.k_printer(parallel_or_singllel, CPU_CORES, E);

    // -------- Part 3:  Solve the problem --------
    // construct solver and solve the case
    let mut eqs: LinearEqs<{ R * C * F }> = LinearEqs::new(
        part.nodes_displacement(),
        part.nodes_force(),
        zero_disp_index,
        part.k(parallel_or_singllel, CPU_CORES).clone(),
    );

    // 1) solve the linear equations of static system using direct method.
    let _ = eqs.solve(calc_method, calc_accuracy, CPU_CORES);

    let calc_time: std::time::Duration = eqs.solver_time_consuming.unwrap();

    // write the displacement and force result into Node2D's field
    part.write_result(&eqs);

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

    // -------- Part 5:  Write clac result into txt file --------
    let problem_type = "BarFrame";
    let element_type = "Rod2D2N";
    let output_path = "/home/zhm/Documents/Scripts/Rust/zhmfem/results/";
    let output_txt = format!(
        "{output_path}{problem_type}_{element_type}_{calc_method}_{parallel_or_singllel}.txt"
    );
    let output_vtk = format!(
        "{output_path}{problem_type}_{element_type}_{calc_method}_{parallel_or_singllel}.vtk"
    );

    part.txt_writer(
        &output_txt,
        calc_time,
        E,
        (strain_energy, external_force_work, potential_energy),
    )
    .expect(">>> !!! Failed to output text result file !!!");

    // Output Calculation result into vtk file
    part.vtk_writer(&output_vtk, element_type)
        .expect(">>> !!! Failed to output vtk file!");

    let total_time = time_start.elapsed();
    println!("\n\n>>> Total time consuming: {:?}", total_time);
}
