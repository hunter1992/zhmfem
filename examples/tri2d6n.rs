use std::collections::HashMap;
use std::time::Instant;

use zhmfem::*;

fn main() {
    // set timing start
    let time_start = Instant::now();

    // set material parameters
    let thick = 1.0 as Dtype;
    let material = (1.0 as Dtype, 0.25 as Dtype); //Young's modulus & Poisson's ratio

    const W: Dtype = 1.0; // width
    const H: Dtype = 1.0; // height

    // number of nodes and freedom
    const R: usize = 3; // rows of nodes
    const C: usize = 3; // columns of nodes
    const M: usize = 6; // node num in single element
    const F: usize = 2; // freedom num in single node

    // construct the solid and mesh it
    //let solid1 = plane::Rectangle::new([0.0 as Dtype, 0.0 as Dtype], [W, H]);
    //let (coords, cpld) = solid1.mesh_with_tri2d6n(R - 1, C - 1);
    let coords: Vec<Vec<Dtype>> = vec![
        vec![1., 0.],
        vec![0., 0.],
        vec![1., 1.],
        vec![0., 1.],
        vec![0.5, 0.],
        vec![0.5, 0.5],
        vec![0., 0.5],
        vec![1., 0.5],
        vec![0.5, 1.],
    ];
    let cpld: Vec<Vec<usize>> = vec![vec![1, 0, 3, 4, 5, 6], vec![3, 0, 2, 5, 7, 8]];

    // set boundary conditions and loads
    let zero_disp: Vec<usize> = vec![2, 3, 6];
    let force_index: Vec<usize> = vec![0, 4];
    let force_value: Vec<Dtype> = vec![-1.0, 1.0];
    let force_data: HashMap<usize, Dtype> = force_index
        .into_iter()
        .zip(force_value.into_iter())
        .collect();
    let nodes: Vec<Node2D> = nodes2d_vec(&coords, &force_data, false);

    let mut tri2d6n_vec: Vec<Tri2D6N> = tri2d6n_vec(thick, &nodes, &cpld);

    let mut part1: Part2D<Tri2D6N, { R * C }, F, M> =
        Part2D::new(1, &nodes, &mut tri2d6n_vec, &cpld, &material);

    let mut eqs: LinearEqs<{ R * C * F }> =
        LinearEqs::new(part1.disps(), part1.forces(), zero_disp, *part1.k(material));

    part1.k_printer(0.);
    eqs.lu_direct_solver();
    part1.write_result(&eqs);

    print_1darr("qe", &part1.disps(), 0.0);
    print_1darr("fe", &part1.forces(), 0.0);

    part1.elems[0].print_strain([1., 0., 0.]);
    part1.elems[0].print_stress([1., 0., 0.], material);

    println!("\n>>> System energy:");
    println!("\tE_d: {:-9.6} (deform energy)", part1.strain_energy());
    println!("\tW_f: {:-9.6} (Exforce works)", part1.force_work());
    println!(
        "\tE_p: {:-9.6} (potential energy)",
        part1.potential_energy()
    );
}
