use std::collections::HashMap;
use zhmfem::*;

fn main() {
    // ------ Part 1: set parameters ------
    // set material parameters
    let material = (1.0f64, 0.25f64); // Young's modulus and Poisson's ratio

    // rectangular simulation area parameters
    let thick = 1.0f64;

    // number of nodes and freedom
    const R: usize = 2; // rows of nodes
    const C: usize = 2; // columns of nodes
    const M: usize = 4; // node num in single element
    const F: usize = 2; // freedom num in single node

    // ------ Part 2: set problems ------
    // construct the solid and mesh it
    let coords = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![1.0, 1.0],
        vec![0.0, 1.0],
    ];
    let cpld = vec![vec![1, 2, 3, 4]];

    // set boundary conditions and mesh it
    // The following two lines are set for
    // a specific problem,not general code
    let zero_disp: Vec<usize> = vec![0, 1, 6];
    let force_index: Vec<usize> = vec![2, 4];
    let force_value: Vec<f64> = vec![-1.0, 1.0];
    let force_data: HashMap<usize, f64> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // tramsform points into nodes
    let nodes: Vec<Node2D> = nodes2d_vec(&coords, &zero_disp, &force_data);

    // construct elements by coupled nodes
    let mut rects: Vec<Quad2D4N> = quad2d4n_vec(thick, &nodes, &cpld);

    // assemble global stiffness matrix
    let mut p1: Part2D<Quad2D4N, { R * C }, F, M> = Part2D::new(1, &nodes, &mut rects, &cpld);
    p1.k(material);
    p1.k_printer(0.0);

    // construct solver and solve the case
    let mut solver: Solver<{ R * C * F }> = Solver::new(p1.disps(), p1.forces(), *p1.k(material));
    solver.solve_static_lu();

    p1.write_result(&solver);

    print_1darr("qe", &p1.disps());
    print_1darr("fe", &p1.forces());

    println!(">>> System energy:");
    println!("\tE_d: {:-9.6} (deform energy)", p1.strain_energy());
    println!("\tW_f: {:-9.6} (exforce works)", p1.force_work());
    println!("\tE_p: {:-9.6} (potential energy)", p1.potential_energy());

    for i in rects.iter() {
        println!("{}", i);
        i.k_printer(0.0);
    }

    println!();
    let target1 = [1.0, -1.0];
    print_1darr("Disp at (1, 0)", &rects[0].point_disp(target1));
    print_1darr("Strain at (1, 0)", &rects[0].strain(target1));
    print_1darr("Stress at (1, 0)", &rects[0].stress(target1, material));

    let target2 = [0.8, 0.8];
    print_1darr("Disp at (0.8, 0.8)", &rects[0].point_disp(target2));
    print_1darr("Strain at (0.8, 0.8)", &rects[0].strain(target2));
    print_1darr("Stress at (0.8,0.8)", &rects[0].stress(target2, material));

    let target3 = [0.8, 0.8];
    print_1darr("Disp at (0.8, 0.8)", &rects[0].point_disp(target3));
    print_1darr("Strain at (0.8, 0.8)", &rects[0].strain(target3));
    print_1darr("Stress at (0.8,0.8)", &rects[0].stress(target3, material));
}
