extern crate nalgebra as na;

use na::Vector3;
use zhmfem::*;

fn main() {
    run();
}

fn run() {
    let ee = 1.0f64;
    let nu = 0.25f64;
    let t = 1.0f64;
    let parameters = (ee, nu, t);

    let node1 = Node2D::new(1, [0.0, 0.0]);
    let node2 = Node2D::new(2, [1.0, 0.0]);
    let node3 = Node2D::new(3, [0.0, 1.0]);
    let node4 = Node2D::new(4, [1.0, 1.0]);

    let tri1 = Triangle::new(1, [&node1, &node2, &node3]);
    tri1.info();
    let tri1_k = tri1.calc_k(parameters);
    println!("{:?}", tri1_k);

    let rect1 = Rectangle::new(2, [&node1, &node2, &node3, &node4]);

    // try to calculate tri1's stiffness matrix K
    //let tri1_k = tri1.k(parameters);

    let vec1 = Vector3::new(1.0, 2.0, 3.0);
    println!("{}", vec1);
}
