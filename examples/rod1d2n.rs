use std::collections::HashMap;
use std::io::{BufWriter, Write};
use zhmfem::*;

fn main() {
    let section_area = 1.0f64;
    let material = (8.0f64, 0.25f64);

    let points = vec![vec![0.0], vec![1.0]];
    let cpld = vec![vec![1, 2]];
    let zero_disp: Vec<usize> = vec![0];
    let force_idx: Vec<usize> = vec![1];
    let force_vlu: Vec<f64> = vec![1.0];
    let force_data: HashMap<usize, f64> =
        force_idx.into_iter().zip(force_vlu.into_iter()).collect();

    let nodes: Vec<Node1D> = nodes1d_vec(&points, &zero_disp, &force_data);

    let mut rod_vec: Vec<Rod1D2N> = rod1d2n_vec(section_area, &nodes, &cpld);

    rod_vec[0].k(material);
    rod_vec[0].k_printer(0.0);
    print_1darr("strain:", &rod_vec[0].strain());

    /*
    let filename = "/home/zhm/Desktop/test_rod1d2n.txt";
    let file = std::fs::File::create(filename).unwrap();
    let mut writer = BufWriter::new(file);
    write!(writer, "{}", rod1.info()).expect("Write failed!");
    */
}
