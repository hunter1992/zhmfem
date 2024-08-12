fn main() {
    println!("");
    println!("      ========================================");
    println!("      ========== WELCOME to zhmfem! ==========");
    println!("      ========================================\n");
    println!(">>> check the examples under 'zhmfem/examples/' path using:");
    println!("        cargo run --examples <example-name>");

    let a: [[zhmfem::Dtype; 2]; 2] = [[1., 2.], [3., 4.]];
    let b: [[zhmfem::Dtype; 3]; 2] = [[10., 20., 30.], [40., 50., 60.]];
    let c: [[zhmfem::Dtype; 2]; 3] = [[100., 200.], [300., 400.], [500., 600.]];

    let z1 = zhmfem::matrix_hstack(&a, &b);
    let z2 = zhmfem::matrix_vstack(&a, &c);
    zhmfem::print_2darr("z1", &z1, 0.);
    zhmfem::print_2darr("z2", &z2, 0.);
}
