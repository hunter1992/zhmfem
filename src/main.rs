use std::collections::HashMap;
use zhmfem::*;

fn main() {
    run();
}

fn run() {
    // set material parameters
    let thick = 1.0f64;
    let material = (1.0f64, 0.25f64);

    // initial input data
    const R: usize = 5;
    const C: usize = 5;
    const M: usize = 3;
    const F: usize = 2;
    let r1 = plane::Rectangle::new([0.0, 0.0], [1.0, 1.0]);
    let (coords, cpld) = r1.mesh_with_tri(R, C);

    let zero_disp: Vec<usize> = vec![0, 1, 40];
    let force_index: Vec<usize> = vec![8, 48];
    let force_value: Vec<f64> = vec![-1.0, 1.0];
    let force_data: HashMap<usize, f64> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // transform points into nodes
    let nodes: Vec<Node2D> = nodes2d_vec(&coords, &zero_disp, &force_data);

    // construct element by coupled nodes
    let mut tris: Vec<Tri2D3N> = tri2d3n_vec(thick, &nodes, &cpld);
    //for i in tris.iter_mut() {
    //    println!("{}", i);
    //    i.k_printer(material);
    //}

    // assemble global stiffness matrix
    let mut p1: Part2D<Tri2D3N, { R * C }, F, M> = Part2D::new(1, &nodes, &mut tris, &cpld);
    //print_2darr("K", p1.k(material));

    //print_1darr("qe", &p1.disps());
    //print_1darr("fe", &p1.forces());

    // construct solver and solve the result
    let mut solver: Solver<{ R * C * F }> = Solver::new(p1.disps(), p1.forces(), *p1.k(material));
    solver.solve_static();

    p1.write_result(&solver);

    //print_1darr("qe", &p1.disps());
    //print_1darr("fe", &p1.forces());

    println!("Part[{}] deform energy: {}", p1.id, p1.strain_energy());
    println!("Part[{}] force work: {}", p1.id, p1.force_work());
    println!(
        "Part[{}] potential energy: {}",
        p1.id,
        p1.potential_energy()
    );
}
