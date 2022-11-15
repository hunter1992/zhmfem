use crate::K;

// two generic const:
// N for N_NODE, F for N_FREEDOM
pub struct Part<Elem: K, const N: usize, const F: usize>
where
    [[f64; N * F]; N * F]: Sized,
{
    pub id: usize,
    pub elems: Vec<Elem>,
    pub cplds: Vec<Vec<usize>>,
    pub k_matrix: Option<[[f64; N * F]; N * F]>,
}

impl<Elem: K, const N: usize, const F: usize> Part<Elem, N, F>
where
    [[f64; N * F]; N * F]: Sized,
{
    pub fn new(id: usize, elems: Vec<Elem>, cplds: Vec<Vec<usize>>) -> Part<Elem, N, F>
    where
        [[f64; N * F]; N * F]: Sized,
    {
        Part {
            id,
            elems,
            cplds,
            k_matrix: None,
        }
    }

    pub fn calc_k(&mut self, material: (f64, f64, f64)) -> [[f64; N * F]; N * F] {
        let mut kk = [[0.0; N * F]; N * F];

        if self.cplds.len() != self.elems.len() {
            println!("---> Error! From assemble_global_k func.");
            println!("     The count of elements not equal to K mat size.");
            panic!("---> Global K failed!");
        }

        // 计算并缓存每个单元的刚度矩阵
        let ks: Vec<_> = self.elems.iter_mut().map(|x| x.k(material)).collect();

        // 用单个单元的节点数，遍历单元的局部刚度矩阵
        let n_nodes: Vec<usize> = (0..self.cplds[0].len()).collect();
        let loc: Vec<Vec<usize>> = full_combination(&n_nodes);

        // 整体刚度矩阵中需要修改的节点坐标对 的位置
        // 注意！这种写法默认cplds字段中节点编号从1开始
        let loc_g: Vec<Vec<Vec<usize>>> = self.cplds.iter().map(|x| full_combination(&x)).collect();

        for i in 0..loc_g.len() {
            for j in 0..loc.len() {
                for k in 0..F {
                    for l in 0..F {
                        kk[(loc_g[i][j][0] - 1) * F + k][(loc_g[i][j][1] - 1) * F + l] += (self.elems[i].k(material))[loc[j][0] * F + k][loc[j][1] * F + l];
                    }
                }
            }
        }

        kk
    }
}

// 生成aim中元素的所有组合
pub fn full_combination(aim: &Vec<usize>) -> Vec<Vec<usize>> {
    let mut rlt: Vec<Vec<usize>> = Vec::new();
    for i in 0..aim.len() {
        for j in 0..aim.len() {
            rlt.push(vec![aim[i], aim[j]]);
        }
    }
    rlt
}
