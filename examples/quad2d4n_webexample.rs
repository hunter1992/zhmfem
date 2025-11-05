#![allow(incomplete_features)]
#![feature(generic_const_exprs)]

use std::thread;
use std::time::Instant;

use zhmfem::*;

fn main() {
    const STACK_SIZE: usize = 512 * 1024 * 1024; // MB = 1024 * 1024
    let executor = thread::Builder::new()
        .stack_size(STACK_SIZE)
        .spawn(|| {
            // set timing start
            let time_start = Instant::now();

            // -------- Part 0: Set initial parameters --------
            const _E: Dtype = 0.0; // Exponent in scientific notation to base 10
            const CPU_CORES: usize = 2;

            // "lu"       for LU       decomposition algorithm or
            // "cholesky" for Cholesky decomposition algorithm or
            // "gs"       for gauss-seidel iteration algorithm
            let calc_method: &str = "cholesky";
            // Calculation accuracy of iterative algorithm
            let calc_accuracy: Dtype = 0.001;

            // "s" or "singllel" or "p" or "parallel"
            let parallel_or_singllel: &str = "s";

            let thick: Dtype = 10.0; //Thickness of the plate
            let material: [Dtype; 2] = [210000.0, 0.30]; //Young's modulud & Poisson's ratio

            // -------- Part 1:  Meshing and applying boundary conditions --------
            // Set mesh and freedom parameters
            const R: usize = 41; // rows of nodes
            const C: usize = 41; // columns of nodes
            const D: usize = R * C; // Dimension of part stiffness matrix
            const M: usize = 4; // num of nodes in single element
            const F: usize = 2; // num of degree freedom at single node

            // Automatically set coords and grouped nodes index
            // Auto-mesh generate coords and grouped nodes index
            const W: Dtype = 1000.0; // width
            const H: Dtype = 1000.0; // height
            let solid1 = Rectangle::new([0.0 as Dtype, 0.0 as Dtype], [W, H]);
            let (points, grpdnidx) = solid1.mesh_with_quad2d4n(R, C);

            // Set boundary conditions and external loads automatically
            let zero_disp_index: Vec<usize> = vec![0, 1, C * (R - 1) * F];
            let force_index: Vec<usize> = vec![(C - 1) * F, (C * R - 1) * F];
            let force_value: Vec<Dtype> = vec![-10000000.0, 10000000.0];

            // -------- Part 2:  Construct nodes, elements and parts --------
            // Construct 2D nodes vector
            let nodes = nodes2d_vec(&points, &force_index, &force_value);

            // Construct Tri2D3N elements vector
            let mut quads = quad2d4n_vec(thick, &nodes, &grpdnidx, material);

            // Construct 2D part & assembly global stiffness matrix
            let mut part: Part2D<'_, Quad2D4N<'_>, D, F, M> =
                Part2D::new(1, &nodes, &mut quads, &grpdnidx);
            // part.k_printer(parallel_or_singllel, CPU_CORES, E);

            // -------- Part 3:  Solve the problem --------
            // construct solver and solve the case
            let mut eqs: LinearEqs<{ R * C * F }> = LinearEqs::new(
                part.nodes_displacement(),
                part.nodes_force(),
                zero_disp_index,
                part.k(parallel_or_singllel, CPU_CORES).clone(),
            );

            // 1) solve the linear equations of static system using direct method.
            // eqs.lu_direct_solver(); //LU decomposition method
            // let output_file = "LU.txt";

            // 2) solve the linear equations of static system using iter method.
            // eqs.gauss_seidel_iter_solver(0.001);
            // let output_file = "G-S.txt";

            // 3) or you can solve the problem with a more concise call:
            eqs.solve(calc_method, calc_accuracy);

            let _calc_time: std::time::Duration = eqs.solver_time_consuming.unwrap();

            // write the displacement and force result into Node2D's field
            part.write_result(&eqs);

            // -------- Part 4:  Print all kinds of result --------
            // print_1darr("qe", &part.nodes_displacement(), E, "v");
            // print_1darr("fe", &part.nodes_force(), E, "v");

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
            let problem_type = "WEBexample";
            let element_type = "Quad2D4N";
            let output_path = "/home/zhm/Documents/Scripts/Rust/zhmfem/results/";
            let output_vtk = format!(
        "{output_path}{problem_type}_{element_type}_{calc_method}_{parallel_or_singllel}.vtk"
    );

            // Output Calculation result into vtk file
            part.vtk_writer(&output_vtk, element_type)
                .expect(">>> !!! Failed to output vtk file!");

            let total_time = time_start.elapsed();
            println!("\n>>> Total time consuming: {:?}", total_time);
        })
        .unwrap();

    executor.join().unwrap();
}
