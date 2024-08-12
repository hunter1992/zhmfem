use std::collections::HashMap;
use std::time::Instant;
use zhmfem::*;

fn main() {
    // set time start
    let time_start = Instant::now();

    // ------ Part 1: set parameters ------
    // set material parameters
    let material = (1.0 as Dtype, 0.25 as Dtype); // Young's modulus and Poisson's ratio

    // rectangular simulation area parameters
    let thick = 1.0 as Dtype;

    const W: Dtype = 1.0;
    const H: Dtype = 1.0;

    // number of nodes and freedom
    const R: usize = 2; // rows of nodes
    const C: usize = 2; // columns of nodes
    const M: usize = 4; // node num in single element
    const F: usize = 2; // freedom num in single node

    // ------ Part 2: set problems ------
    // construct the solid and mesh it
    let solid1 = plane::Rectangle::new([0.0 as Dtype, 0.0 as Dtype], [W, H]);
    let (coords, cpld) = solid1.mesh_with_rect(R, C);
    print_2dvec("coords", &coords, 0.);
    print!("cpld{:?}", &cpld);

    // set boundary conditions and mesh it
    // The following two lines are set for
    // a specific problem,not general code
    let zero_disp: Vec<usize> = vec![0, 1, 4];
    let force_index: Vec<usize> = vec![2, 6];
    let force_value: Vec<Dtype> = vec![-1.0, 1.0];
    let force_data: HashMap<usize, Dtype> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // tramsform points into nodes
    let nodes: Vec<Node2D> = nodes2d_vec(&coords, &force_data, false);

    // construct elements by coupled nodes
    let mut rects: Vec<Quad2D4N> = quad2d4n_vec(thick, &nodes, &cpld);

    // assemble global stiffness matrix
    let mut part1: Part2D<Quad2D4N, { R * C }, F, M> =
        Part2D::new(1, &nodes, &mut rects, &cpld, &material);
    part1.k(material);
    part1.k_printer(0.0);

    // construct solver and solve the case
    let mut eqs: LinearEqs<{ R * C * F }> =
        LinearEqs::new(part1.disps(), part1.forces(), zero_disp, *part1.k(material));
    eqs.lu_direct_solver();
    //eqs.gauss_seidel_iter_solver(0.0001);

    part1.write_result(&eqs);

    print_1darr("qe", &part1.disps(), 0.0);
    print_1darr("fe", &part1.forces(), 0.0);

    let stress_point: [[Dtype; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0], [1.0, 1.0]];
    for pt in stress_point.iter() {
        print!("\npt[{}, {}]:", pt[0], pt[1]);
        part1.elems[0].print_strain(*pt);
        part1.elems[0].print_stress(*pt, material);
    }

    println!("\n>>> System energy:");
    println!("\tE_d: {:-9.6} (deform energy)", part1.strain_energy());
    println!("\tW_f: {:-9.6} (exforce works)", part1.force_work());
    println!(
        "\tE_p: {:-9.6} (potential energy)",
        part1.potential_energy()
    );

    let total_time = time_start.elapsed();
    println!("\n>>> Total time consuming: {:?}", total_time);
}
