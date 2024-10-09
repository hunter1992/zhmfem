extern crate nalgebra as na;

use crate::{calc::LinearEqs, node::Node1D, Dtype, K};
use na::*;

/// Three generic const: N for N_NODE, F for N_FREEDOM, M for N_NODE in 1 element
pub struct Part1D<'a, Elem: K, const N: usize, const F: usize, const M: usize>
where
    [[Dtype; N * F]; N * F]: Sized,
{
    pub id: usize,
    pub nodes: &'a mut [Node1D],
    pub elems: &'a mut [Elem],
    pub cplds: &'a [Vec<usize>],
    pub material: &'a (Dtype, Dtype),
    pub k_matrix: Option<[[Dtype; N * F]; N * F]>,
}

impl<'a, Elem: K + 'static, const N: usize, const F: usize, const M: usize>
    Part1D<'a, Elem, N, F, M>
where
    [[Dtype; N * F]; N * F]: Sized,
{
    pub fn new(
        id: usize,
        nodes: &'a mut [Node1D],
        elems: &'a mut [Elem],
        cplds: &'a [Vec<usize>],
        material: &'a (Dtype, Dtype),
    ) -> Self
    where
        [[Dtype; N * F]; N * F]: Sized,
    {
        Part1D {
            id,
            nodes,
            elems,
            cplds,
            material,
            k_matrix: None,
        }
    }

    /// Get displacement of all nodes
    pub fn disps(&self) -> [Dtype; N * F] {
        let mut data: [Dtype; N * F] = [0.0; N * F];
        for idx in 0..N {
            data[idx] = self.nodes[idx].displs[0];
        }
        data
    }

    /// Get force of all nodes
    pub fn forces(&self) -> [Dtype; N * F] {
        let mut data: [Dtype; N * F] = [0.0; N * F];
        for idx in 0..N {
            data[idx] = self.nodes[idx].forces[0];
        }
        data
    }

    /// Get the deform energy of the part
    pub fn strain_energy(&self) -> Dtype {
        if self.k_matrix.is_none() {
            panic!("---> Error! Part[{}]'s stiffness matrix is null!", self.id);
        }
        let q = SMatrix::<Dtype, { N * F }, 1>::from(self.disps());
        let k = SMatrix::<Dtype, { N * F }, { N * F }>::from(self.k_matrix.unwrap()).transpose();
        let rlt: [[Dtype; 1]; 1] = (q.transpose() * k * q).into();
        0.5 * rlt[0][0]
    }

    /// Get the external force work on the part
    pub fn force_work(&self) -> Dtype {
        let q = SMatrix::<Dtype, { N * F }, 1>::from(self.disps());
        let f = SMatrix::<Dtype, { N * F }, 1>::from(self.forces());
        let rlt: [[Dtype; 1]; 1] = (f.transpose() * q).into();
        rlt[0][0]
    }

    pub fn potential_energy(&self) -> Dtype {
        self.strain_energy() - self.force_work()
    }

    /// Write the disp and force result into nodes
    pub fn write_result(&mut self, slv: &LinearEqs<{ N * F }>) {
        let disp = slv.disps;
        let force = slv.forces;
        for (idx, node) in self.nodes.iter_mut().enumerate() {
            node.displs[0] = disp[idx];
            node.forces[0] = force[idx];
        }
    }

    /// Assemble the global stiffness mat K
    pub fn k(&mut self, material: (Dtype, Dtype)) -> &[[Dtype; N * F]; N * F]
    where
        <Elem as K>::Kmatrix: std::ops::Index<usize, Output = [Dtype; M * F]>,
    {
        if self.k_matrix.is_none() {
            if self.cplds.len() != self.elems.len() {
                println!("---> Error! From Part2D.calc_k func.");
                println!("     The count of elements not eq to K mat size.");
                panic!("---> Assembly global K failed!");
            }

            println!(
                "\n>>> Assembling Part1D#{}'s global stiffness matrix K{} ......",
                self.id, self.id
            );
            let mut part_k: [[Dtype; N * F]; N * F] = [[0.0; N * F]; N * F];

            // 计算并缓存每个单元的刚度矩阵
            let elem_ks: Vec<_> = self.elems.iter_mut().map(|x| x.k(material)).collect();

            // 获取单个单元内的节点数目n，构造0到n的range
            // 用于遍历单个单元的局部刚度矩阵k
            let n_nodes_in_elem: Vec<usize> = (0..self.cplds[0].len()).collect();
            let loc: Vec<Vec<usize>> = full_combination(&n_nodes_in_elem);

            // 构造整体刚度矩阵中需要修改的节点坐标对
            // 注意！这种写法默认传进来的coupled_nodes中节点编号从0起
            let loc_g: Vec<Vec<Vec<usize>>> =
                self.cplds.iter().map(|x| full_combination(&x)).collect();

            for i in 0..loc_g.len() {
                for j in 0..loc.len() {
                    for k in 0..F {
                        for l in 0..F {
                            part_k[(loc_g[i][j][0]) * F + k][(loc_g[i][j][1]) * F + l] +=
                                elem_ks[i][loc[j][0] * F + k][loc[j][1] * F + l];
                        }
                    }
                }
            }
            self.k_matrix.get_or_insert(part_k)
        } else {
            self.k_matrix.as_ref().unwrap()
        }
    }

    /// Print part's global stiffness matrix
    pub fn k_printer(&mut self, n_exp: Dtype) {
        if self.k_matrix.is_none() {
            println!(
                "!!! Part #{}'s K matrix is empty! call k() to calc it.",
                self.id
            );
        }

        print!("Part #{}  K =  (* 10^{})\n[", self.id, n_exp as u8);
        for row in 0..N * F {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..N * F {
                print!(
                    " {:>-10.6} ",
                    self.k_matrix.unwrap()[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                );
            }
            if row == (N * F - 1) {
                println!("]]");
            } else {
                println!("]");
            }
        }
        println!("");
    }
}

fn full_combination(aim: &Vec<usize>) -> Vec<Vec<usize>> {
    let mut rlt: Vec<Vec<usize>> = Vec::new();
    for i in 0..aim.len() {
        for j in 0..aim.len() {
            rlt.push(vec![aim[i], aim[j]]);
        }
    }
    rlt
}
