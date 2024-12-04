use crate::{Dtype, Export, LinearEqs, Node2D, K};
use std::io::{BufWriter, Write};

pub struct Part2D<'part2d, Elem: K, const N: usize, const F: usize, const M: usize>
where
    [[Dtype; N * F]; N * F]: Sized,
{
    pub id: usize,
    pub nodes: &'part2d [Node2D],
    pub elems: &'part2d mut [Elem],
    pub cplds: &'part2d [Vec<usize>], //cplds: coupled nodes index
    pub k_matrix: Option<[[Dtype; N * F]; N * F]>,
}

impl<'part2d, Elem: K, const N: usize, const F: usize, const M: usize>
    Part2D<'part2d, Elem, N, F, M>
where
    [[Dtype; N * F]; N * F]: Sized,
{
    pub fn new(
        id: usize,
        nodes: &'part2d Vec<Node2D>,
        elems: &'part2d mut Vec<Elem>,
        cplds: &'part2d Vec<Vec<usize>>,
    ) -> Self {
        Part2D {
            id,
            nodes,
            elems,
            cplds,
            k_matrix: None,
        }
    }

    /// Get displacement of all nodes
    pub fn nodes_displacement(&self) -> [Dtype; N * F] {
        let mut data: [Dtype; N * F] = [0.0; N * F];
        for idx in 0..N {
            data[idx * 2] = self.nodes[idx].displs.borrow()[0];
            data[idx * 2 + 1] = self.nodes[idx].displs.borrow()[1];
        }
        data
    }

    /// Get force of all nodes
    pub fn nodes_force(&self) -> [Dtype; N * F] {
        let mut data: [Dtype; N * F] = [0.0; N * F];
        for idx in 0..N {
            data[idx * 2] = self.nodes[idx].forces.borrow()[0];
            data[idx * 2 + 1] = self.nodes[idx].forces.borrow()[1];
        }
        data
    }

    pub fn k(&mut self) -> &[[Dtype; N * F]; N * F]
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
                "\n>>> Assembling Part2D(#{})'s global stiffness matrix K{} ......",
                self.id, self.id
            );
            let mut part_k: [[Dtype; N * F]; N * F] = [[0.0; N * F]; N * F];

            // 计算并缓存每个单元的刚度矩阵
            let elem_ks: Vec<_> = self.elems.iter_mut().map(|x| x.k()).collect();

            // 获取单个单元内的节点数目n，构造0到n的range，用于遍历单个单元的局部刚度矩阵k
            // 暂时认为part中只有一种单元类型，所有单元内节点数目相同，用cplds[0].len
            let n_nodes_in_elem: Vec<usize> = (0..self.cplds[0].len()).collect();
            let node_loc: Vec<Vec<usize>> = full_combination(&n_nodes_in_elem);

            // 构造整体刚度矩阵中需要修改的节点坐标对
            // 注意！这种写法默认传进来的coupled_nodes中节点编号从1起
            let loc_g: Vec<Vec<Vec<usize>>> =
                self.cplds.iter().map(|x| full_combination(&x)).collect();

            for i in 0..loc_g.len() {
                for j in 0..node_loc.len() {
                    for k in 0..F {
                        for l in 0..F {
                            part_k[(loc_g[i][j][0]) * F + k][(loc_g[i][j][1]) * F + l] +=
                                elem_ks[i][node_loc[j][0] * F + k][node_loc[j][1] * F + l];
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
    pub fn k_printer(&mut self, n_exp: Dtype)
    where
        <Elem as K>::Kmatrix: std::ops::Index<usize, Output = [Dtype; M * F]>,
    {
        if self.k_matrix.is_none() {
            self.k();
        }

        print!("\nPart #{}  K =  (* 10^{})\n[", self.id, n_exp as u8);
        for row in 0..(N * F) {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..(N * F) {
                print!(
                    " {:>-12.6} ",
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

    /// Write the disp and force result into nodes
    pub fn write_result(&mut self, slv: &LinearEqs<{ N * F }>) {
        let disp = slv.disps;
        let force = slv.forces;
        for (idx, node) in self.nodes.iter().enumerate() {
            node.displs.borrow_mut()[0] = disp[idx * 2];
            node.displs.borrow_mut()[1] = disp[idx * 2 + 1];
            node.forces.borrow_mut()[0] = force[idx * 2];
            node.forces.borrow_mut()[1] = force[idx * 2 + 1];
        }
    }
}

impl<'a, Elem: K, const N: usize, const F: usize, const M: usize> Export
    for Part2D<'a, Elem, N, F, M>
where
    [[Dtype; N * F]; N * F]: Sized,
{
    fn txt_writer(
        &self,
        target_file: &str,
        calc_time: std::time::Duration,
        n_exp: Dtype,
        energy: (Dtype, Dtype, Dtype),
    ) -> std::io::Result<bool> {
        let txt_file = std::fs::File::create(target_file).unwrap();
        let mut text_writer = BufWriter::new(txt_file);

        println!("\n>>> Writing calc results into txt file ......");
        write!(text_writer, ">>> ZHMFEM calculating results:").expect("Write txt file error!");

        let (deform, exwork, potential) = energy;
        write!(text_writer, "\n").expect("Write parts' result into yxy file failed!!!");
        write!(text_writer, "\n>>> System energy:\n")
            .expect("Write parts' result into yxy file failed!!!");
        write!(text_writer, "\t\t\tE_d: {:-12.6} (deform energy)\n", deform)
            .expect("Write parts' result into yxy file failed!!!");
        write!(text_writer, "\t\t\tW_f: {:-12.6} (exforce works)\n", exwork)
            .expect("Write parts' result into yxy file failed!!!");
        write!(
            text_writer,
            "\t\t\tE_p: {:-12.6} (potential energy)\n",
            potential
        )
        .expect("Write parts' result into yxy file failed!!!");

        write!(
            text_writer,
            "\n>>> Solver time consuming: {:?}\n",
            calc_time
        )
        .expect("Write txt file error!");

        write!(text_writer, "\n>>> Details of each element:")
            .expect("Write parts' result into txt file failed!!!");
        for elem in self.elems.iter() {
            write!(text_writer, "{}\n", elem.info(n_exp)).expect("Write info failed!");
        }
        text_writer.flush().expect("!!! Flush txt file failed!");
        println!("    Down!");
        Ok(true)
    }

    fn vtk_writer(&self, _target_file: &str) -> std::io::Result<bool> {
        Ok(true)
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
