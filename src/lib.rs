#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(test)]

extern crate nalgebra as na;
extern crate test;

mod calc;
mod elem;
mod mesh;
mod node;
mod part;

pub use atomic_float::{AtomicF32, AtomicF64};
pub use calc::LinearEqs;
pub use elem::dim1::{beam::Beam1D2N, rod::Rod1D2N};
pub use elem::dim2::{
    quadrila::Quad2D4N,
    rod::Rod2D2N,
    triangle::{Tri2D3N, Tri2D6N},
};
pub use mesh::plane::Rectangle;
pub use na::SMatrix;
pub use node::{Node1D, Node2D, Node3D};
pub use part::{part1d::Part1D, part2d::Part2D};

pub use std::collections::HashMap;
pub use std::io::{BufWriter, Write};

pub type Dtype = f32;
pub type ADtype = AtomicF32;
//pub type Dtype = f64;
//pub type ADtype = AtomicF64;

pub type Jacobian2D = SMatrix<Dtype, 2, 2>;

#[derive(Copy, Clone)]
pub enum SS {
    Dim1([Dtype; 1]),
    Dim2([Dtype; 3]),
    Dim3([Dtype; 6]),
}

/// K trait generate element's stiffness matrix under linear analysis.
/// Output stress/strain vector at some point in element using Vector,
/// Kmatrix is the element's stiffness matrix to get.
pub trait K {
    type Kmatrix;

    fn k(&mut self) -> &Self::Kmatrix
    where
        Self::Kmatrix: std::ops::Index<usize>;
    fn k_printer(&self, n_exp: Dtype);
    fn k_string(&self, n_exp: Dtype) -> String;

    fn strain_intpt(&mut self) -> &SS;
    fn stress_intpt(&mut self) -> &SS;

    fn info(&self, n_exp: Dtype) -> String;
    fn id(&self) -> usize;
}

/// Export trait generate Part's informations and write it into files
/// Export trait generate txt or vtk files for output
pub trait Export {
    fn txt_writer(
        &self,
        target_file: &str,
        calc_time: std::time::Duration,
        n_exp: Dtype,
        energy: (Dtype, Dtype, Dtype),
    ) -> std::io::Result<bool>;

    fn vtk_writer(&mut self, target_file: &str, elem_type: &str) -> std::io::Result<bool>;
}

/// Formatted print 1d array with scientific form
pub fn print_1darr<const C: usize>(name: &str, arr: &[Dtype; C], n_exp: Dtype, h_or_v: &str) {
    println!("\n{} = (10^{} *)", name, n_exp);
    match h_or_v {
        "h" => {
            print!("[[");
            for col in 0..C {
                if col == 0 {}
                print!(
                    " {:-18.6} ",
                    arr[col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                );
            }
            println!("]]\n");
        }
        "v" => {
            print!("[[\n");
            for row in 0..C {
                if row == 0 {}
                print!(
                    " {:-18.6} \n",
                    arr[row] / (10.0_f64.powf(n_exp as f64)) as Dtype
                );
            }
            println!("]]\n");
        }
        _ => panic!("!!! Wrong print arg for print_1darr, from /src/lib/fn print_1darr !!!"),
    }
}

/// Formated print a 2D array with scientific form
pub fn print_2darr<const R: usize, const C: usize>(
    name: &str,
    id: usize,
    array: &[[Dtype; C]; R],
    n_exp: Dtype,
) {
    println!("\n{}[{}] = (10^{} *)", name, id, n_exp);
    for row in 0..R {
        if row == 0 {
            print!("[[");
        } else {
            print!(" [");
        }
        for col in 0..C {
            print!(
                " {:>-18.6} ",
                array[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
            );
        }
        if row == array.len() - 1 {
            println!("]]\n");
        } else {
            println!("]");
        }
    }
}

/// Formated print a 1D vector with scientific form
pub fn print_1dvec(name: &str, vec: &[Dtype], n_exp: Dtype, h_or_v: &str) {
    println!("\n{} = (10^{} *)", name, n_exp);
    print!("[[");
    if h_or_v == "h" {
        for &ele in vec.iter() {
            print!(" {:-18.6} ", ele / (10.0_f64.powf(n_exp as f64)) as Dtype);
        }
    } else if h_or_v == "v" {
        for &ele in vec.iter() {
            print!(" {:-18.6} \n", ele / (10.0_f64.powf(n_exp as f64)) as Dtype);
        }
    }
    println!("]]\n");
}

/// Formated print a 2D vector with scientific form
pub fn print_2dvec(name: &str, mat: &[Vec<Dtype>], n_exp: Dtype) {
    println!("\n{} = (10^{} *)", name, n_exp);
    for row in 0..mat.len() {
        if row == 0 {
            print!("[[");
        } else {
            print!(" [");
        }
        for col in 0..mat[0].len() {
            print!(
                " {:-18.6} ",
                mat[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
            );
        }
        if row == mat.len() - 1 {
            println!("]]\n");
        } else {
            println!("]");
        }
    }
}

/// Constructe a 1D nodes vector
pub fn nodes1d_vec(coords: &[Vec<Dtype>], forces: &HashMap<usize, Dtype>) -> Vec<Node1D> {
    let mut nodes: Vec<Node1D> = Vec::with_capacity(coords.len());
    for (idx, coord) in coords.iter().enumerate() {
        nodes.push(Node1D::new(idx, [coord[0]]));
    }
    for (idx, &f) in forces {
        nodes[idx / 1].forces.borrow_mut()[idx % 1] = f;
    }
    nodes
}

/// convert 1D or 2D points to 3D points
pub fn convert_to_3d_points<'zhmfem>(
    dim: usize,
    points: &'zhmfem mut Vec<Vec<Dtype>>,
) -> &'zhmfem Vec<Vec<Dtype>> {
    let l: usize = points.len();
    if dim == 1 {
        for idx in 0..l {
            points[idx].push(0.0);
            points[idx].push(0.0);
        }
    } else if dim == 2 {
        for idx in 0..l {
            points[idx].push(0.0);
        }
    } else {
        panic!("!!! convert 1d/2d points to 3d points failed.");
    }
    points
}

/// Constructe a 2D nodes vector
pub fn nodes2d_vec(coords: &[Vec<Dtype>], forces: &HashMap<usize, Dtype>) -> Vec<Node2D> {
    let mut nodes: Vec<Node2D> = Vec::with_capacity(coords.len());
    for (idx, coord) in coords.iter().enumerate() {
        nodes.push(Node2D::new(idx, [coord[0], coord[1]]));
    }
    for (idx, &f) in forces {
        nodes[idx / 2].forces.borrow_mut()[idx % 2] = f;
    }
    nodes
}

/// Constructe a 3D nodes vector
pub fn nodes3d_vec(coords: &[Vec<Dtype>], forces: &HashMap<usize, Dtype>) -> Vec<Node3D> {
    let mut nodes: Vec<Node3D> = Vec::with_capacity(coords.len());
    for (idx, coord) in coords.iter().enumerate() {
        nodes.push(Node3D::new(idx, [coord[0], coord[1], coord[2]]));
    }
    for (idx, &f) in forces {
        nodes[idx / 3].forces.borrow_mut()[idx % 3] = f;
    }
    nodes
}

/// calculate part's deform energy
pub fn strain_energy<const D: usize>(
    stiffness_matrix: [[Dtype; D]; D],
    displacement: [Dtype; D],
) -> Dtype {
    let disp = SMatrix::<Dtype, D, 1>::from(displacement);
    let k_matrix = SMatrix::<Dtype, D, D>::from(stiffness_matrix);
    let strain_energy: [[Dtype; 1]; 1] = (0.5 * disp.transpose() * k_matrix * disp).into();
    strain_energy[0][0]
}

/// calculate external forces' work
pub fn external_force_work<const D: usize>(
    external_force: [Dtype; D],
    displacement: [Dtype; D],
) -> Dtype {
    let disp = SMatrix::<Dtype, D, 1>::from(displacement);
    let external_force = SMatrix::<Dtype, D, 1>::from(external_force);
    let strain_energy: [[Dtype; 1]; 1] = (external_force.transpose() * disp).into();
    strain_energy[0][0]
}

/// calculate part's potential energy
pub fn potential_energy<const D: usize>(
    stiffness_matrix: [[Dtype; D]; D],
    external_force: [Dtype; D],
    displacement: [Dtype; D],
) -> Dtype {
    strain_energy(stiffness_matrix, displacement)
        - external_force_work(external_force, displacement)
}

/// Constructe a rod1d2n elements vector
/// every rod with same cross sectional area
pub fn rod1d2n_vec<'rod1d2n>(
    //nodes: &'rod1d2n Vec<Node1D>,
    nodes: &'rod1d2n [Node1D],
    coupled_nodes_idx: &[Vec<usize>],
    cross_sectional_area: &'rod1d2n [Dtype],
    material: &'rod1d2n (Dtype, Dtype),
) -> Vec<Rod1D2N<'rod1d2n>> {
    let mut rod1d2n: Vec<Rod1D2N> = Vec::with_capacity(coupled_nodes_idx.len());
    for (ele_idx, cpld) in coupled_nodes_idx.iter().enumerate() {
        rod1d2n.push(Rod1D2N::new(
            ele_idx,
            cross_sectional_area[ele_idx],
            [&nodes[cpld[0]], &nodes[cpld[1]]],
            material,
        ))
    }
    rod1d2n
}

/// Constructe a rod2d2n elements vector
/// every rod with same cross sectional area
pub fn rod2d2n_vec<'rod2d2n>(
    nodes: &'rod2d2n [Node2D],
    coupled_nodes_idx: &[Vec<usize>],
    cross_sectional_area: &'rod2d2n [Dtype],
    material: &'rod2d2n (Dtype, Dtype),
) -> Vec<Rod2D2N<'rod2d2n>> {
    let mut rod2d2n: Vec<Rod2D2N> = Vec::with_capacity(coupled_nodes_idx.len());
    for (ele_idx, cpld) in coupled_nodes_idx.iter().enumerate() {
        rod2d2n.push(Rod2D2N::new(
            ele_idx,
            cross_sectional_area[ele_idx],
            [&nodes[cpld[0]], &nodes[cpld[1]]],
            material,
        ))
    }
    rod2d2n
}

/// Constructe a beam1d2n elements vector
/// every beam with same cross sectional area
pub fn beam1d2n_vec<'beam1d2n>(
    nodes: &'beam1d2n [Node2D],
    coupled_nodes_idx: &[Vec<usize>],
    moment_of_inertia: &'beam1d2n [Dtype],
    cross_sectional_area: &'beam1d2n [Dtype],
    material: &'beam1d2n (Dtype, Dtype),
) -> Vec<Beam1D2N<'beam1d2n>> {
    let mut beam1d2n: Vec<Beam1D2N> = Vec::with_capacity(coupled_nodes_idx.len());
    for (ele_idx, cpld) in coupled_nodes_idx.iter().enumerate() {
        beam1d2n.push(Beam1D2N::new(
            ele_idx,
            moment_of_inertia[ele_idx],
            cross_sectional_area[ele_idx],
            [&nodes[cpld[0]], &nodes[cpld[1]]],
            material,
        ))
    }
    beam1d2n
}

/// Constructe a tri2d3n elements vector
/// every element with same thick
pub fn tri2d3n_vec<'tri2d3n>(
    thick: Dtype,
    nodes: &'tri2d3n [Node2D],
    coupled_nodes_idx: &[Vec<usize>],
    material: &'tri2d3n (Dtype, Dtype),
) -> Vec<Tri2D3N<'tri2d3n>> {
    let mut tri2d3n: Vec<Tri2D3N> = Vec::with_capacity(coupled_nodes_idx.len());
    for (ele_idx, cpld) in coupled_nodes_idx.iter().enumerate() {
        tri2d3n.push(Tri2D3N::new(
            ele_idx,
            thick,
            [&nodes[cpld[0]], &nodes[cpld[1]], &nodes[cpld[2]]],
            material,
        ))
    }
    tri2d3n
}

/// Constructe a quadrila elements vector
/// every element with same thick
pub fn quad2d4n_vec<'quad2d4n>(
    thick: Dtype,
    nodes: &'quad2d4n Vec<Node2D>,
    coupled_nodes: &[Vec<usize>],
    material: &'quad2d4n (Dtype, Dtype),
) -> Vec<Quad2D4N<'quad2d4n>> {
    let mut quad2d4n: Vec<Quad2D4N> = Vec::with_capacity(coupled_nodes.len());
    for (ele_idx, cpld) in coupled_nodes.iter().enumerate() {
        quad2d4n.push(Quad2D4N::new(
            ele_idx,
            thick,
            [
                &nodes[cpld[0]],
                &nodes[cpld[1]],
                &nodes[cpld[2]],
                &nodes[cpld[3]],
            ],
            material,
        ))
    }
    quad2d4n
}

#[cfg(test)]
mod testing {
    use super::*;
    use test::Bencher;

    #[test]
    fn gen_nodes() {
        let node1 = Node1D::new(1, [1.0]);
        let node2 = Node2D::new(2, [1.0, 2.0]);
        let node3 = Node3D::new(3, [1.0, 2.0, 3.0]);

        assert_eq!(1usize, node1.id);
        assert_ne!(3usize, node2.id);
        assert_eq!(3usize, node3.id);

        assert_eq!([1.0 as Dtype], node1.coords);
        assert_eq!([1.0 as Dtype, 2.0 as Dtype], node2.coords);
        assert_eq!([1.0 as Dtype, 2.0 as Dtype, 3.0 as Dtype], node3.coords);
    }

    #[test]
    fn gen_elements() {
        let material: (Dtype, Dtype) = (1.0, 0.25);
        // 1D
        let node_a = Node1D::new(1, [0.0]);
        let node_b = Node1D::new(2, [1.0]);
        let node_c = Node1D::new(3, [3.0]);

        let rod1 = Rod1D2N::new(1, 1.0, [&node_a, &node_b], &material);
        let rod2 = Rod1D2N::new(2, 1.0, [&node_b, &node_c], &material);
        assert_eq!(1usize, rod1.id);
        assert_ne!(2usize, rod1.id);
        assert_eq!(2usize, rod2.id);

        // 2D
        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [0.0, 1.0]);
        let node3 = Node2D::new(3, [1.0, 0.0]);
        let node4 = Node2D::new(4, [1.0, 1.0]);
        let thick = 1.0;

        let tri1 = Tri2D3N::new(1, thick, [&node1, &node2, &node3], &material);
        let tri2 = Tri2D3N::new(2, thick, [&node4, &node2, &node3], &material);

        let rec1 = Quad2D4N::new(3, thick, [&node1, &node2, &node3, &node4], &material);

        assert_eq!(1usize, tri1.id);
        assert_ne!(2usize, tri1.id);
        assert_eq!(3usize, rec1.id);

        assert_ne!([10.0 as Dtype, 11.0 as Dtype], tri2.nodes[0].coords);
        assert_eq!([1.0 as Dtype, 1.0 as Dtype], tri2.nodes[0].coords);

        assert_eq!(vec![0.0, 0.0, 1.0], tri1.get_nodes_xcoords());
        assert_eq!(vec![0.0, 0.0, 1.0, 1.0], rec1.get_nodes_xcoords());
        assert_ne!(vec![0.0, 0.0, 1.0], tri1.get_nodes_ycoords());
        assert_ne!(vec![0.0, 0.0, 1.0, 1.0], rec1.get_nodes_ycoords());

        assert_eq!(0.5 as Dtype, tri1.area());
        assert_eq!(0.5 as Dtype, tri2.area());
        assert_eq!(1.0 as Dtype, rec1.area());
    }

    #[test]
    fn gen_parts() {
        let material = (1.0 as Dtype, 0.25 as Dtype);
        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [1.0, 0.0]);
        let node3 = Node2D::new(3, [1.0, 1.0]);
        let node4 = Node2D::new(4, [0.0, 1.0]);
        let thick = 1.0;

        let nodes = vec![node1, node2, node3, node4];
        let cplds = vec![vec![0, 1, 3], vec![1, 2, 3]];
        let mut tris: Vec<Tri2D3N> = tri2d3n_vec(thick, &nodes, &cplds, &material);

        let p1: Part2D<Tri2D3N, 4, 2, 3> = Part2D::new(1, &nodes, &mut tris, &cplds);
        assert_eq!(p1.elems[1].nodes[1].coords[1], 1.0);
        assert_ne!(p1.elems[1].nodes[1].coords[1], -1.0);

        let disp = vec![
            -1024., -1024., -1024., -1024., -1024., -1024., -1024., -1024.,
        ];
        let force = vec![0., 0., 0., 0., 0., 0., 0., 0.];
        let nodes_disp = p1.nodes_displacement();
        let nodes_force = p1.nodes_force();
        assert_ne!(disp, nodes_disp);
        assert_eq!(force, nodes_force);
    }

    #[test]
    fn mesh_rect_with_tri() {
        // set rect's width and height
        const W: Dtype = 1.0;
        const H: Dtype = 1.0;

        // number of nodes and freedom
        const R: usize = 2; //rows of nodes
        const C: usize = 2; //rows of nodes

        let rect_geo = Rectangle::new([0., 0.], [W, H]);
        let (coords, coupled_nodes) = rect_geo.mesh_with_tri2d3n(R, C);

        assert_eq!(coords, [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]);
        assert_eq!(coupled_nodes, [[0, 1, 2], [3, 2, 1]]);
        /*   2____3
         *   |\   |
         *   | \  |
         *   |  \ |
         *   |   \|
         *   0____1
         */
    }

    #[test]
    fn calc_tri_elem_k() {
        let material = (1.0 as Dtype, 0.25 as Dtype);

        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [1.0, 0.0]);
        let node3 = Node2D::new(3, [1.0, 1.0]);
        let node4 = Node2D::new(4, [0.0, 1.0]);
        let thick = 1.0;

        let mut tri1 = Tri2D3N::new(1, thick, [&node1, &node2, &node4], &material);
        let mut tri2 = Tri2D3N::new(2, thick, [&node3, &node4, &node2], &material);

        let k1 = tri1.k();
        let k2 = tri2.k();

        assert_eq!(k1, k2);
    }

    #[test]
    fn calc_quad_elem_k() {
        let material = (1.0 as Dtype, 0.25 as Dtype);

        let thick = 1.0;
        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [1.0, 0.0]);
        let node3 = Node2D::new(3, [1.0, 1.0]);
        let node4 = Node2D::new(4, [0.0, 1.0]);

        let mut quad1 = Quad2D4N::new(1, thick, [&node1, &node2, &node3, &node4], &material);

        let k1 = quad1.k();
        assert_eq!(0.48888892 as Dtype, k1[0][0]);
    }

    #[bench]
    /// benchmark的结果是:393.49 ns +/- 36.23 ns/iter (Intel 8265U 插电) 20241205
    fn calc_quad_k_speed(b: &mut Bencher) {
        b.iter(|| calc_quad_elem_k());
    }

    #[bench]
    /// benchmark的结果是:229.13 ns +/- 06.28 ns/iter (Intel 8265U 插电) 20241205
    fn calc_tri_k_speed(b: &mut Bencher) {
        b.iter(|| calc_tri_elem_k());
    }
}
