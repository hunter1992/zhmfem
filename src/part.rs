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

    pub fn cache_elem_ks(&mut self, material: (f64, f64, f64)){
        if self.cplds.len() != self.elems.len() {
            println!("---> Error! From assemble_global_k func.");
            println!("     The count of elements not equal to K mat size.");
            panic!("---> Global K failed!");
        }
        // 计算并缓存每个单元的刚度矩阵
        let _ks: Vec<_> = self.elems.iter_mut().map(|x| x.k(material)).collect();
    }

    pub fn assemble_global_k(&mut self, material:(f64,f64,f64))->
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
