#![feature(generic_const_exprs)]

mod elem;
mod node;

pub use elem::*;
pub use node::*;

pub fn full_combination(aim: &Vec<usize>) -> Vec<Vec<usize>> {
    let mut rlt: Vec<Vec<usize>> = Vec::new();
    for i in 0..aim.len() {
        for j in 0..aim.len() {
            rlt.push(vec![aim[i], aim[j]]);
        }
    }
    rlt
}

pub fn print_vec2d(mat: &[Vec<f64>]) {
    print!("[");
    for row in 0..mat.len() {
        if row == 0 {
            print!("[");
        } else {
            print!(" [");
        }
        for col in 0..mat[0].len() {
            print!(" {:-9.6} ", mat[row][col]);
        }
        if row == mat.len() - 1 {
            println!("]]");
        } else {
            println!("]");
        }
    }
}

pub fn print_arr2d<T, const R: usize, const C: usize>(arr: &[[T; C]; R])
where
    T: std::fmt::Display,
{
    for r in 0..R {
        if r == 0 {
            print!("[[");
        } else {
            print!(" [");
        }
        for c in 0..C {
            print!(" {:-9.6} ", arr[r][c]);
        }
        if r == arr.len() - 1 {
            println!("]]");
        } else {
            println!("]");
        }
    }
}

pub fn nodes1d_vec(points: &[Vec<f64>]) -> Vec<Node1D> {
    let mut nodes: Vec<Node1D> = Vec::new();
    for (idx, point) in points.iter().enumerate() {
        nodes.push(Node1D::new(idx + 1, [point[0]]));
    }
    nodes
}

pub fn nodes2d_vec(points: &[Vec<f64>]) -> Vec<Node2D> {
    let mut nodes: Vec<Node2D> = Vec::new();
    for (idx, point) in points.iter().enumerate() {
        nodes.push(Node2D::new(idx + 1, [point[0], point[1]]));
    }
    nodes
}

pub fn nodes3d_vec(points: &[Vec<f64>]) -> Vec<Node3D> {
    let mut nodes: Vec<Node3D> = Vec::new();
    for (idx, point) in points.iter().enumerate() {
        nodes.push(Node3D::new(idx + 1, [point[0], point[1], point[2]]));
    }
    nodes
}

pub fn tri2d3n_vec<'a>(nodes: &'a [Node2D], couples: &'a [Vec<usize>]) -> Vec<Triangle<'a>> {
    let mut tri2d3n: Vec<Triangle> = Vec::new();
    for (ele_id, cpld) in couples.iter().enumerate() {
        tri2d3n.push(Triangle::new(
            ele_id + 1,
            [
                &nodes[cpld[0] - 1],
                &nodes[cpld[1] - 1],
                &nodes[cpld[2] - 1],
            ],
        ));
    }
    tri2d3n
}

pub fn global_k<const N_NODES: usize, const N_FREEDOM: usize>(
    material: (f64, f64, f64),
    coupled_nodes: &Vec<Vec<usize>>,
    elem_vec: &mut Vec<Triangle>,
) -> [[f64; N_NODES * N_FREEDOM]; N_FREEDOM * N_NODES] {
    let mut kk = [[0.0; N_NODES * N_FREEDOM]; N_FREEDOM * N_NODES];

    if coupled_nodes.len() != elem_vec.len() {
        println!("---> Error! From assemble_global_k func.");
        println!("     The count of elements not equal to K mat size.");
        panic!("---> Global K failed!");
    }

    // 计算并缓存每个单元的刚度矩阵
    let ks: Vec<_> = elem_vec.iter_mut().map(|x| x.k(material)).collect();

    // 获取单个单元内的节点数目N,构造0到N的range
    // 用于遍历单个单元的局部刚度矩阵k
    let n_nodes: Vec<usize> = (0..coupled_nodes[0].len()).collect();
    let loc: Vec<Vec<usize>> = full_combination(&n_nodes);

    // 整体刚度矩阵中需要修改的节点坐标对
    // 注意！这种写法默认传进来的coupled_nodes中节点编号从1起
    let loc_g: Vec<Vec<Vec<usize>>> = coupled_nodes.iter().map(|x| full_combination(&x)).collect();

    for i in 0..loc_g.len() {
        for j in 0..loc.len() {
            for k in 0..N_FREEDOM {
                for l in 0..N_FREEDOM {
                    kk[(loc_g[i][j][0] - 1) * N_FREEDOM + k]
                        [(loc_g[i][j][1] - 1) * N_FREEDOM + l] = kk
                        [(loc_g[i][j][0] - 1) * N_FREEDOM + k]
                        [(loc_g[i][j][1] - 1) * N_FREEDOM + l]
                        + ks[i][loc[j][0] * N_FREEDOM + k][loc[j][1] * N_FREEDOM + l];
                }
            }
        }
    }
    kk
}

#[cfg(test)]
mod testing {
    use super::*;

    #[test]
    fn gen_nodes() {
        let node1 = Node1D::new(1, [1.0]);
        let node2 = Node2D::new(2, [1.0, 2.0]);
        let node3 = Node3D::new(3, [1.0, 2.0, 3.0]);

        assert_eq!(1usize, node1.id);
        assert_ne!(3usize, node2.id);
        assert_eq!(3usize, node3.id);

        assert_eq!([1.0f64, 2.0f64], node2.coord);
        assert_eq!([1.0f64, 2.0f64, 3.0f64], node3.coord);
        assert_eq!([1.0f64], node1.coord);
    }

    #[test]
    fn gen_elements() {
        let node1 = Node2D::new(1, [0.0, 0.0]);
        let node2 = Node2D::new(2, [0.0, 1.0]);
        let node3 = Node2D::new(3, [1.0, 0.0]);
        let node4 = Node2D::new(4, [1.0, 1.0]);

        let tri1 = Triangle::new(1, [&node1, &node2, &node3]);
        let tri2 = Triangle::new(2, [&node4, &node2, &node3]);

        let rec1 = Rectangle::new(3, [&node1, &node2, &node3, &node4]);

        assert_eq!(1usize, tri1.id);
        assert_ne!(2usize, tri1.id);
        assert_eq!(3usize, rec1.id);

        assert_ne!([10.0f64, 11.0f64], tri2.nodes[0].coord);
        assert_eq!([1.0f64, 1.0f64], tri2.nodes[0].coord);

        assert_eq!(vec![0.0, 0.0, 1.0], tri1.xs());
        assert_eq!(vec![0.0, 0.0, 1.0, 1.0], rec1.get_xs());
        assert_ne!(vec![0.0, 0.0, 1.0], tri1.ys());
        assert_ne!(vec![0.0, 0.0, 1.0, 1.0], rec1.get_ys());

        assert_eq!(0.5f64, tri1.area());
        assert_eq!(0.5f64, tri2.area());
    }
}
