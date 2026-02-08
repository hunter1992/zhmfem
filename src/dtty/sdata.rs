use crate::dtty::basic::Dtype;
use std::fmt;

#[derive(Clone)]
pub struct NodeSData2D {
    pub sdata_on_node: Vec<[Dtype; 3]>,
    pub area_arround_node: Vec<Dtype>,
}

impl Default for NodeSData2D {
    fn default() -> Self {
        NodeSData2D {
            sdata_on_node: vec![[0., 0., 0.]],
            area_arround_node: vec![0.0],
        }
    }
}

impl NodeSData2D {
    pub fn new(sdata_on_node: Vec<[Dtype; 3]>, area_arround_node: Vec<Dtype>) -> Self {
        NodeSData2D {
            sdata_on_node,
            area_arround_node,
        }
    }

    pub fn data_smoothing(&self, average_method: &str) -> [Dtype; 3] {
        assert_eq!(self.sdata_on_node.len(), self.area_arround_node.len());
        match average_method {
            "a" => self.arithmetic_mean(),
            "arithmetic" => self.arithmetic_mean(),
            "w" => self.weighted_mean(),
            "weighted" => self.weighted_mean(),
            "auto" => self.weighted_mean(),
            _ => panic!("Please provide an algorithm for strain/stress smoothing."),
        }
    }

    fn arithmetic_mean(&self) -> [Dtype; 3] {
        assert_eq!(self.sdata_on_node.len(), self.area_arround_node.len());

        let mut s11: Dtype = 0.0;
        let mut s22: Dtype = 0.0;
        let mut s12: Dtype = 0.0;
        self.sdata_on_node
            .iter()
            .zip(self.area_arround_node.iter())
            .map(|(sdata_vec, _)| {
                s11 += sdata_vec[0];
                s22 += sdata_vec[1];
                s12 += sdata_vec[2];
            })
            .count();
        let inv: Dtype = 1.0 / ((self.area_arround_node.len() - 1) as Dtype);
        [s11 * inv, s22 * inv, s12 * inv]
    }

    fn weighted_mean(&self) -> [Dtype; 3] {
        assert_eq!(self.sdata_on_node.len(), self.area_arround_node.len());

        let mut s11: Dtype = 0.0;
        let mut s22: Dtype = 0.0;
        let mut s12: Dtype = 0.0;
        self.sdata_on_node
            .iter()
            .zip(self.area_arround_node.iter())
            .map(|(sdata_vec, area)| {
                s11 += sdata_vec[0] * area;
                s22 += sdata_vec[1] * area;
                s12 += sdata_vec[2] * area;
            })
            .count();
        let inv: Dtype = self.area_arround_node.iter().sum::<Dtype>().recip();
        [s11 * inv, s22 * inv, s12 * inv]
    }
}

impl fmt::Display for NodeSData2D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            ">>> NodeSdata:\n\tsdatas_on_node:\n\t{:?}\n\tarea_arround_node:\n\t{:?}\n",
            self.sdata_on_node, self.area_arround_node
        )
    }
}
