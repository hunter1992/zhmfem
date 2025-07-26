use crate::dtty::basic::Dtype;

/// Store node strain or stress data
#[derive(Clone)]
pub struct Sdata {
    pub sdata: Vec<Vec<Dtype>>,
    pub id_info: Vec<usize>,   //这个点被哪几个单元共有
    pub area_info: Vec<Dtype>, //共有这个点的几个单元的面积
}

impl Sdata {
    pub fn new(sdata: Vec<Vec<Dtype>>, id_info: Vec<usize>, area_info: Vec<Dtype>) -> Self {
        Sdata {
            sdata,
            id_info,
            area_info,
        }
    }
}
