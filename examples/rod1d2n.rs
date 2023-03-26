use std::collections::HashMap;
use std::io::{BufWriter, Write};
use zhmfem::*;

fn main() {
    let section_area = 1.0f64;
    let material = (8.0f64, 0.25f64);

    let node1 = Node1D::new(1, [0.0]);
    let node2 = Node1D::new(2, [1.0]);

    let points = vec![vec![0.0], vec![1.0]];
    let cpld = vec![vec![1, 2]];
    let zero_disp: Vec<usize> = vec![0];
    let force_idx: Vec<usize> = vec![1];
    let force_vlu: Vec<f64> = vec![1.0];
    let force_data: HashMap<usize, f64> =
        force_idx.into_iter().zip(force_vlu.into_iter()).collect();

    let nodes: Vec<Node1D> = nodes1d_vec(&points, &zero_disp, &force_data);

    for i in nodes {
        println!("{}", i);
    }

    /*
    let filename = "/home/zhm/Desktop/test_rod1d2n.txt";
    let file = std::fs::File::create(filename).unwrap();
    let mut writer = BufWriter::new(file);
    write!(writer, "{}", rod1.info()).expect("Write failed!");
    */
}
