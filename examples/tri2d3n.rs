use std::collections::HashMap;
use zhmfem::*;

fn main() {
    // set material parameters
    let material = (1.0f64, 0.25f64); //Young's modulud & Poisson's ratio

    // rectangular simulation area parameters
    let thick = 1.0f64;
    const W: f64 = 1.0; // width
    const H: f64 = 1.0; // height

    // number of nodes and freedom
    const R: usize = 2; // rows of nodes
    const C: usize = 2; // columns of nodes
    const M: usize = 3; // node num in single element
    const F: usize = 2; // freedom num in single node

    // construct the solid and mesh it
    let solid1 = plane::Rectangle::new([0.0, 0.0], [W, H]);
    let (coords, cpld) = solid1.mesh_with_tri(R, C);

    // set boundary conditions and loads
    let zero_disp: Vec<usize> = vec![0, 1, (R - 1) * C * 2];
    let force_index: Vec<usize> = vec![(C - 1) * 2, (R * C - 1) * 2];
    let force_value: Vec<f64> = vec![-1.0, 1.0];
    let force_data: HashMap<usize, f64> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();

    // transform points into nodes
    let nodes: Vec<Node2D> = nodes2d_vec(&coords, &zero_disp, &force_data);

    // construct elements by coupled nodes
    let mut tris: Vec<Tri2D3N> = tri2d3n_vec(thick, &nodes, &cpld);

    // assemble global stiffness matrix
    let mut p1: Part2D<Tri2D3N, { R * C }, F, M> = Part2D::new(1, &nodes, &mut tris, &cpld);
    print_2darr("\nK", p1.k(material));

    // construct solver and solve the case
    let mut solver: Solver<{ R * C * F }> = Solver::new(p1.disps(), p1.forces(), *p1.k(material));
    solver.solve_static();

    p1.write_result(&solver);

    print_1darr("qe", &p1.disps());
    print_1darr("fe", &p1.forces());

    println!("deform energy: {}", p1.strain_energy());
    println!("force work: {}", p1.force_work());
    println!("potential energy: {}", p1.potential_energy());

    for i in tris.iter() {
        println!("{}", i);
        i.print_strain();
        i.print_stress(material);
    }
    print_1darr("\nDisp at (1, 0)", &tris[0].point_disp([1.0, 0.0]));
    print_1darr("\nDisp at (1, 1)", &tris[1].point_disp([1.0, 1.0]));
}
