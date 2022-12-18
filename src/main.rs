use std::collections::HashMap;
use zhmfem::*;

fn main() {
    println!("");
    println!("      ========================================");
    println!("      ========== WELCOME TO zhmfem! ==========");
    println!("      ========================================\n");
    println!(">>> check the examples under 'zhmfem/examples/' path using:");
    println!("        cargo run --examples <example-name>");

    let thick = 1.0f64;
    let material = (1.0f64, 0.25f64);
    //let thick = 0.1f64;
    //let material = (30000000.0f64, 0.3f64);

    //let coords: Vec<Vec<f64>> = vec![vec![1., 0.], vec![2., 0.], vec![2.25, 1.5], vec![1.25, 1.]];
    let coords: Vec<Vec<f64>> = vec![vec![0., 0.], vec![1., 0.], vec![1.0, 1.0], vec![0.0, 1.0]];
    let cplds: Vec<Vec<usize>> = vec![vec![1, 2, 3, 4]];

    let zero_disp: Vec<usize> = vec![0, 1, 6];
    let force_idx: Vec<usize> = vec![2, 4];
    let force_vle: Vec<f64> = vec![-1., 1.];
    let force_data: HashMap<usize, f64> =
        force_idx.into_iter().zip(force_vle.into_iter()).collect();

    let nodes: Vec<Node2D> = nodes2d_vec(&coords, &zero_disp, &force_data);
    let mut rects: Vec<Quad2D4N> = rec2d4n_vec(thick, &nodes, &cplds);

    println!("{}", rects[0]);
    print_2darr("k1", rects[0].k(material));

    let mut p1: Part2D<Quad2D4N, 4, 2, 4> = Part2D::new(1, &nodes, &mut rects, &cplds);
    let mut solver: Solver<8> = Solver::new(p1.disps(), p1.forces(), *p1.k(material));
    solver.solve_static();
    p1.write_result(&solver);

    print_1darr("qe", &p1.disps());
    print_1darr("fe", &p1.forces());

    println!("deform energy: {}", p1.strain_energy());
    println!("force work: {}", p1.force_work());
    println!("potential energy: {}", p1.potential_energy());
}
