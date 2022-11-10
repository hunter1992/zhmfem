use crate::elem::K;

// two generic const:
// N for N_NODE, F for N_FREEDOM
pub struct Part<Elem: K, const N: usize, const F: usize>
where
    [[f64; N * F]; N * F]: Sized,
{
    pub id: usize,
    pub elems: Vec<Elem>,
    pub coupled_nodes: Vec<Vec<usize>>,
    pub k_matrix: Option<[[f64; N * F]; N * F]>,
}

impl<Elem: K, const N: usize, const F: usize> Part<Elem, N, F>
where
    [[f64; N * F]; N * F]: Sized,
{
    pub fn new(id: usize, elems: Vec<Elem>, coupled_nodes: Vec<Vec<usize>>) -> Part<Elem, N, F>
    where
        [[f64; N * F]; N * F]: Sized,
    {
        Part {
            id,
            elems,
            coupled_nodes,
            k_matrix: None,
        }
    }
}

// 生成aim中元素的所有组合
//pub fn full_combination(aim: &Vec<usize>) -> Vec<Vec<usize>> {
//    let mut rlt: Vec<Vec<usize>> = Vec::new();
//    for i in 0..aim.len() {
//        for j in 0..aim.len() {
//            rlt.push(vec![aim[i], aim[j]]);
//        }
//    }
//    rlt
//}
