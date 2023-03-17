use std::io::{BufWriter, Write};
use zhmfem::*;

fn main() {
    let section_area = 1.0f64;
    let material = (8.0f64, 0.25f64);

    let node1 = Node1D::new(1, [0.0]);
    let node2 = Node1D::new(2, [1.0]);

    let mut rod1 = Rod1D2N::new(1, section_area, [&node1, &node2]);

    println!("{}", rod1);
    rod1.k(material);
    rod1.k_printer(0.0);

    let filename = "/home/zhm/Desktop/test_rod1d2n.txt";
    let file = std::fs::File::create(filename).unwrap();
    let mut writer = BufWriter::new(file);
    write!(writer, "{}", rod1.info()).expect("Write failed!");
}
