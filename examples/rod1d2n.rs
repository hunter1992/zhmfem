// Fundamentals of Finite Element Analysis PAGE17 example2.3(1)
use std::time::Instant;

use zhmfem::*;

fn main() {
    // set time start
    let time_start = Instant::now();

    const E: Dtype = 0.0; // Exponent in scientific notation to base 10
    let cpu_cores: usize = 2;

    // "lu"       for LU       decomposition algorithm or
    // "cholesky" for Cholesky decomposition algorithm or
    // "pardiso"  for calling  Panua Tech's PARDISO
    // "gs"       for gauss-seidel iteration algorithm
    let calc_method: &str = "gs";

    // Calculation accuracy of iterative algorithm
    let calc_accuracy: Dtype = 0.000001;

    // "s" or "singllel" or "p" or "parallel"
    let parallel_or_singllel: &str = "s";

    let section_areas: [Dtype; 3] = [0.02, 0.03, 0.06];
    let material = [200000.0 as Dtype, 0.25 as Dtype];

    const R: usize = 1;
    const C: usize = 4;
    const M: usize = 2;
    const F: usize = 1;

    let points: Vec<[Dtype; 1]> = vec![[0.0], [0.1], [0.2], [0.3]];
    let cpld = vec![vec![0, 1], vec![1, 2], vec![2, 3]];
    let zero_disp_index: Vec<usize> = vec![3];
    let force_idx: Vec<usize> = vec![0, 2];
    let force_vlu: Vec<Dtype> = vec![-100.0, 50.0];

    let nodes: Vec<Node1D> = nodes1d_vec(&points, &force_idx, &force_vlu);

    let mut rod_vec: Vec<Rod1D2N> = rod1d2n_vec(&nodes, &cpld, &section_areas, material);

    let mut part: Part1D<Rod1D2N, { R * C }, F, M> = Part1D::new(1, &nodes, &mut rod_vec, &cpld);

    let mut eqs: LinearEqs<{ R * C * F }> = LinearEqs::new(
        part.nodes_displacement(),
        part.nodes_force(),
        zero_disp_index,
        part.k(parallel_or_singllel, cpu_cores).clone(),
    );

    let _ = eqs.solve(calc_method, calc_accuracy, cpu_cores);

    let calc_time: std::time::Duration = eqs.solver_time_consuming.unwrap();

    part.write_result(&eqs);

    print_1darr("qe", &part.nodes_displacement(), 0.0, "v");
    print_1darr("fe", &part.nodes_force(), 0.0, "v");

    println!("\n>>> System energy:");
    let strain_energy: Dtype = strain_energy(
        part.k(parallel_or_singllel, cpu_cores).clone(),
        part.nodes_displacement(),
    );
    let external_force_work: Dtype =
        external_force_work(part.nodes_force(), part.nodes_displacement());
    let potential_energy: Dtype = potential_energy(
        part.k(parallel_or_singllel, cpu_cores).clone(),
        part.nodes_force(),
        part.nodes_displacement(),
    );
    println!("\tE_d: {:-9.6} (deform energy)", strain_energy);
    println!("\tW_f: {:-9.6} (exforce works)", external_force_work);
    println!("\tE_p: {:-9.6} (potential energy)", potential_energy);

    let problem_type = "one-dim_stretch";
    let element_type = "Rod1D2N";
    let output_path = "/home/zhm/Documents/Scripts/Rust/zhmfem/results/";
    let output_txt = format!(
        "{output_path}{problem_type}_{element_type}_{calc_method}_{parallel_or_singllel}.txt"
    );
    let output_vtk = format!(
        "{output_path}{problem_type}_{element_type}_{calc_method}_{parallel_or_singllel}.vtk"
    );

    // Output Calculation result into txt file
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
