use crate::{full_combination, node::*, K};

/// two generic const:
/// N for N_NODE, F for N_FREEDOM, M for N_NODE in 1 element
pub struct Part2D<Elem: K, const N: usize, const F: usize, const M:usize>
where
    [[f64; N * F]; N * F]: Sized,
{
    pub id: usize,
    pub elems: Vec<Elem>,
    pub cplds: Vec<Vec<usize>>,
    pub material: (f64, f64, f64),
    pub k_matrix: Option<[[f64; N * F]; N * F]>,
}

impl<Elem: K, const N: usize, const F: usize, const M:usize> Part2D<Elem, N, F, M>
where
    [[f64; N * F]; N * F]: Sized,
{
    pub fn new(
        id: usize,
        elems: Vec<Elem>,
        cplds: Vec<Vec<usize>>,
        material: (f64, f64, f64),
    ) -> Self
    where
        [[f64; N * F]; N * F]: Sized,
    {
        Part2D {
            id,
            elems,
            cplds,
            material,
            k_matrix: None,
        }
    }

    /// Get displacement on all nodes
    pub fn get_nodes_disp(&self, nodes: &Vec<Node2D>) -> Vec<f64> {
        let mut data: Vec<f64> = Vec::with_capacity(N * F);
        for node in nodes.iter() {
            data.push(node.disps[0]);
            data.push(node.disps[1]);
        }
        data
    }

    /// Get force on all nodes
    pub fn get_nodes_force(&self, nodes: &Vec<Node2D>) -> Vec<f64> {
        let mut data: Vec<f64> = Vec::with_capacity(N * F);
        for node in nodes.iter() {
            data.push(node.force[0]);
            data.push(node.force[1]);
        }
        data
    }

    /// Assemble the global stiffness mat K
    pub fn k(&mut self) -> &[[f64; N * F]; N * F]
    where
        <Elem as K>::Kmatrix: std::ops::Index<usize, Output = [f64; M * F]>,
    {
        if self.k_matrix.is_none() {
            if self.cplds.len() != self.elems.len() {
                println!("---> Error! From Part2D.calc_k func.");
                println!("     The count of elements not eq to K mat size.");
                panic!("---> Assembly global K failed!");
            }

            let mut part_k: [[f64; N * F]; N * F] = [[0.0; N * F]; N * F];

            // 计算并缓存每个单元的刚度矩阵
            let elem_ks: Vec<_> = self.elems.iter_mut().map(|x| x.k(self.material)).collect();

            // 获取单个单元内的节点数目n，构造0到n的range
            // 用于遍历单个单元的局部刚度矩阵k
            let n_nodes_in_elem: Vec<usize> = (0..self.cplds[0].len()).collect();
            let loc: Vec<Vec<usize>> = full_combination(&n_nodes_in_elem);

            // 构造整体刚度矩阵中需要修改的节点坐标对
            // 注意！这种写法默认传进来的coupled_nodes中节点编号从1起
            let loc_g: Vec<Vec<Vec<usize>>> =
                self.cplds.iter().map(|x| full_combination(&x)).collect();

            for i in 0..loc_g.len() {
                for j in 0..loc.len() {
                    for k in 0..F {
                        for l in 0..F {
                            part_k[(loc_g[i][j][0] - 1) * F + k][(loc_g[i][j][1] - 1) * F + l] +=
                                (elem_ks[i])[loc[j][0] * F + k][loc[j][1] * F + l];
                        }
                    }
                }
            }
            self.k_matrix.get_or_insert(part_k)
        } else {
            self.k_matrix.as_ref().unwrap()
        }
    }
}
