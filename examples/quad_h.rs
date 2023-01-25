use std::collections::HashMap;
use zhmfem::*;

fn main() {
    // ------ Part 1: set parameters ------
    // set material parameters
    let material = (1.0f64, 0.25f64); // Young's modulus and Poisson's ratio

    // quad area
    let thick = 0.1f64;

    // number of nodes and freedom
    const R: usize = 2; // rows of nodes
    const C: usize = 2; // columns of nodes
    const M: usize = 4; // node num in single element
    const F: usize = 2; // freedom num in single node

    // ------ Part 2: set problems ------
    // construct the solid and mesh it
    let points = vec![
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
    let force_index: Vec<usize> = vec![3, 5];
    let force_value: Vec<f64> = vec![1.0, 1.0];
    let force_data: HashMap<usize, f64> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // tramsform points into nodes
    let nodes: Vec<Node2D> = nodes2d_vec(&points, &zero_disp, &force_data);

    // construct elements by coupled nodes
    let mut quads: Vec<Quad2D4N> = quad2d4n_vec(thick, &nodes, &cpld);

    // assemble global stiffness matrix
    let mut p1: Part2D<Quad2D4N, { R * C }, F, M> = Part2D::new(1, &nodes, &mut quads, &cpld);
    p1.k(material);
    p1.k_printer(0.0);

    let mut solver: Solver<{ R * C * F }> = Solver::new(p1.disps(), p1.forces(), *p1.k(material));
    solver.solve_static();
    p1.write_result(&solver);

    print_1darr("qe", &p1.disps());
    print_1darr("fe", &p1.forces());

    for i in quads.iter() {
        println!("{}", i);
        i.k_printer(0.0);
    }

    let target1 = [1.0, -1.0];
    print_1dvec("Disp at (1, 0)", &quads[0].point_disp(target1));
    print_1dvec("Strain at (1, 0)", &quads[0].strain(target1));
    print_1dvec("Stress at (1, 0)", &quads[0].stress(target1, material));

    let target2 = [0.0, 0.0];
    print_1dvec("Disp at (0, 0)", &quads[0].point_disp(target2));
    print_1dvec("Strain at (0, 0)", &quads[0].strain(target2));
    print_1dvec("Stress at (0, 0)", &quads[0].stress(target2, material));
}
