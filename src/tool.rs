use crate::dtty::{basic::Dtype, matrix::CompressedMatrix};
use crate::elem::{
    dim1::{beam::Beam1D2N, rod::Rod1D2N},
    dim2::{quadrila::Quad2D4N, rod::Rod2D2N, triangle::Tri2D3N},
};
use crate::node::{Node1D, Node2D, Node3D};
use na::SMatrix;
use std::collections::HashMap;

/// Return a matrix compressed by Skyline symmetry algorithm
pub fn compress_matrix<const D: usize>(mat: [[Dtype; D]; D]) -> CompressedMatrix {
    let mut value: Vec<Dtype> = vec![];
    let mut ptr: Vec<usize> = vec![];

    let mut flag: bool = true;
    let mut counter: usize = 0;
    for idx in 0..D {
        ptr.push(counter);
        for idy in 0..=idx {
            if mat[idx][idy] == 0.0 {
                if flag == false {
                    // 正定矩阵必满秩，下面处理全零行(或列)的分支暂时注释掉
                    /*if idy == idx {
                        value.push(mat[idx][idy]);
                        counter += 1;
                    }*/
                    continue;
                } else {
                    value.push(mat[idx][idy]);
                    counter += 1;
                    continue;
                }
            }
            value.push(mat[idx][idy]);
            counter += 1;
            flag = true;
        }
        flag = false;
    }
    ptr.push(value.len());

    CompressedMatrix { value, ptr }
}

/// calculate part's deform energy
pub fn strain_energy<const D: usize>(
    stiffness_matrix: CompressedMatrix,
    displacement: [Dtype; D],
) -> Dtype {
    let disp = SMatrix::<Dtype, D, 1>::from(displacement);
    let k_matrix = SMatrix::<Dtype, D, D>::from(stiffness_matrix.recover());
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
    stiffness_matrix: CompressedMatrix,
    external_force: [Dtype; D],
    displacement: [Dtype; D],
) -> Dtype {
    strain_energy(stiffness_matrix, displacement)
        - external_force_work(external_force, displacement)
}

// Formatted print 1d array with scientific form
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
