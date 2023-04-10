use std::collections::HashMap;
use zhmfem::*;

fn main() {
    // set material parameters
    let material = (10000000.0f64, 0.33f64); //Young's modulud & Poisson's ratio

    // rectangular simulation area parameters
    let thick = 0.1f64;

    // number of nodes and freedom
    const R: usize = 2; // rows of nodes
    const C: usize = 2; // columns of nodes
    const M: usize = 3; // node num in single element
    const F: usize = 2; // freedom num in single node

    // construct the solid and mesh it
    let points = vec![vec![2., 1.], vec![2., 0.], vec![0., 1.], vec![0., 0.]];
    let cpld = vec![vec![1, 3, 2], vec![4, 2, 3]];

    // set boundary conditions and loads
    let zero_disp: Vec<usize> = vec![4, 5, 6, 7];
    let force_index: Vec<usize> = vec![1, 3];
    let force_value: Vec<f64> = vec![-50000.0, -50000.0];
    let force_data: HashMap<usize, f64> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // transform points into nodes
    let nodes: Vec<Node2D> = nodes2d_vec(&points, &zero_disp, &force_data);

    // construct elements by coupled nodes
    let mut tris: Vec<Tri2D3N> = tri2d3n_vec(thick, &nodes, &cpld);

    // assemble global stiffness matrix
    let mut p1: Part2D<Tri2D3N, { R * C }, F, M> = Part2D::new(1, &nodes, &mut tris, &cpld);
    println!("");
    p1.k(material);
    p1.k_printer(6.0);

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

    println!("\n==================== ELEMENT INFO ====================");

    for i in tris.iter() {
        println!("{}", i);
        i.k_printer(6.0);
        i.print_strain();
        i.print_stress(material);
    }
}
