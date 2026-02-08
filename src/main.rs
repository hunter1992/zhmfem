fn main() {
    println!("");
    println!("    ========================================");
    println!("    ========== WELCOME to ZHMFEM! ==========");
    println!("    ========================================\n");
    println!(">>> check the examples under 'zhmfem/examples/' path");
    println!("    and run a example  with command:");
    println!();
    println!("        cargo run --example <example-name>");
    println!();
    println!("    or:");
    println!();
    println!("        cargo run -j 4 --release");

    println!("\n\n\n>>> ====== ZHMFEM Data size (Byte)======");
    use std::mem;
    println!("\n    Sizeof Dtype :     {}", size_of::<zhmfem::Dtype>());
    let size_of_node1d = mem::size_of::<zhmfem::Node1D>();
    let size_of_node2d = mem::size_of::<zhmfem::Node2D>();
    let size_of_node3d = mem::size_of::<zhmfem::Node3D>();
    println!(
        "\n    Sizeof Node1D:    {}\n           Node2D:    {}\n           Node3D:    {}",
        size_of_node1d, size_of_node2d, size_of_node3d
    );

    let size_of_rod1d2n = mem::size_of::<zhmfem::Rod1D2N>();
    let size_of_tri2d3n = mem::size_of::<zhmfem::Tri2D3N>();
    let size_of_quad2d4n = mem::size_of::<zhmfem::Quad2D4N>();
    let size_of_rod2d2n = mem::size_of::<zhmfem::Rod2D2N>();
    println!(
        "\n    Sizeof Rod1D2N:   {}\n           Tri2D3N:   {}\n           Quad2D4N:  {}\n           Rod2D2N:   {}",
        size_of_rod1d2n, size_of_tri2d3n, size_of_quad2d4n, size_of_rod2d2n
    );

    let size_of_part2d_rod1d2n = mem::size_of::<zhmfem::Part1D<'_, zhmfem::Rod1D2N<'_>, 2, 1, 2>>();
    let size_of_part2d_tri2d3n = mem::size_of::<zhmfem::Part2D<'_, zhmfem::Tri2D3N<'_>, 4, 2, 3>>();
    let size_of_part2d_quad2d4n =
        mem::size_of::<zhmfem::Part2D<'_, zhmfem::Quad2D4N<'_>, 4, 2, 4>>();
    let size_of_part2d_rod2d2n = mem::size_of::<zhmfem::Part2D<'_, zhmfem::Rod2D2N<'_>, 2, 2, 2>>();
    println!(
        "\n    Sizeof Part1D(Rod1D2N mesh):    {}",
        size_of_part2d_rod1d2n
    );
    println!(
        "           Part2D(Tri2D3N  mesh):   {}",
        size_of_part2d_tri2d3n
    );
    println!(
        "           Part2D(Quad2D4N mesh):   {}",
        size_of_part2d_quad2d4n
    );
    println!(
        "           Part2D(Rod2D2N mesh):    {}",
        size_of_part2d_rod2d2n
    );
    println!("\n=========================================\n");
}
