use std::collections::HashMap;
use std::io::{BufWriter, Write};
use zhmfem::*;

fn main() {
    // set material parameters
    let thick = 1.0f64;
    let material = (1.0f64, 0.25f64); //Young's modulud & Poisson's ratio

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
    let mut tri_vec: Vec<Tri2D3N> = tri2d3n_vec(thick, &nodes, &cpld);

    // assemble global stiffness matrix
    let mut part1: Part2D<Tri2D3N, { R * C }, F, M> = Part2D::new(1, &nodes, &mut tri_vec, &cpld);
    println!("");
    part1.k(material);
    part1.k_printer(0.0);

    // construct solver and solve the case
    let mut solver: Solver<{ R * C * F }> =
        Solver::new(part1.disps(), part1.forces(), *part1.k(material));
    solver.solve_static();

    part1.write_result(&solver);

    print_1darr("qe", &part1.disps());
    print_1darr("fe", &part1.forces());

    println!(">>> System energy:");
    println!("\tE_d: {:-9.6} (deform energy)", part1.strain_energy());
    println!("\tW_f: {:-9.6} (exforce works)", part1.force_work());
    println!(
        "\tE_p: {:-9.6} (potential energy)",
        part1.potential_energy()
    );
    println!("\n==================== ELEMENT INFO ====================");

    for i in tri_vec.iter() {
        println!("{}", i);
        i.k_printer(0.0);
        i.print_strain();
        i.print_stress(material);
    }

    println!();
    print_1darr("Disp at (1, 0)", &tri_vec[0].point_disp([1.0, 0.0]));
    print_1darr("Strain at (1, 0)", &tri_vec[0].strain());
    print_1darr("Stress at (1, 0)", &tri_vec[0].stress(material));

    print_1darr("Disp at (0.8, 0.8)", &tri_vec[0].point_disp([0.8, 0.8]));
    print_1darr("Strain at (0.8, 0.8)", &tri_vec[0].strain());
    print_1darr("Stress at (0.8,0.8)", &tri_vec[0].stress(material));

    print_1darr("Disp at (0.8, 0.8)", &tri_vec[1].point_disp([0.8, 0.8]));
    print_1darr("Strain at (0.8, 0.8)", &tri_vec[1].strain());
    print_1darr("Stress at (0.8,0.8)", &tri_vec[1].stress(material));

    // Write the result into file
    let file_name = "/home/zhm/Desktop/tri2d3n_result.txt";
    let file = std::fs::File::create(file_name).unwrap();
    let mut writer = BufWriter::new(file);
    write!(writer, ">>> ZHMFEM calculating result:").expect("Write error!");

    for elem in tri_vec.iter() {
        write!(writer, "{}", elem.info()).expect("Write info failed!");
        write!(writer, "\tStrain: {:-9.6?}\n", elem.strain()).expect("Write strain failed!");
        write!(writer, "\tStress: {:-9.6?}\n", elem.stress(material))
            .expect("Write stress failed!");
        write!(
            writer,
            "\n\tStiffness matrix k{} = \n{}\n",
            elem.id,
            elem.k_string(0.0) //设置刚度矩阵元素科学记数次数
        )
        .expect("!!! Write k matrix failed!");
    }
    writer.flush().expect("!!! Flush failed!");
}