use std::collections::HashMap;
use zhmfem::*;

fn main() {
    // ------ Part 1: set parameters ------
    // set material parameters
    //let material = (3000000.0f64, 0.3f64); // Young's modulus and Poisson's ratio
    let material = (1.0f64, 0.25f64); // Young's modulus and Poisson's ratio

    // quad area
    let thick = 1.0f64;

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
    let force_index: Vec<usize> = vec![2, 4];
    let force_value: Vec<f64> = vec![-1.0, 1.0];
    let force_data: HashMap<usize, f64> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // transform points into nodes
    let nodes: Vec<Node2D> = nodes2d_vec(&points, &zero_disp, &force_data);

    // construct elements by coupled nodes
    let mut quads: Vec<Quad2D4N> = quad2d4n_vec(thick, &nodes, &cpld);
    print_2darr("J", &quads[0].jacobian([-1., 1.]));

    // assemble global stiffness matrix
    let mut p1: Part2D<Quad2D4N, { R * C }, F, M> = Part2D::new(1, &nodes, &mut quads, &cpld);
    p1.k(material);
    p1.k_printer(0.0);

    let mut eqs: LinearEqs<{ R * C * F }> =
        LinearEqs::new(p1.disps(), p1.forces(), *p1.k(material));
    eqs.lu_direct_solver();
    p1.write_result(&eqs);

    print_1darr("qe", &p1.disps());
    print_1darr("fe", &p1.forces());
}
