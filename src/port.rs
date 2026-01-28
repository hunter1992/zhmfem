use crate::dtty::{basic::Dtype, matrix::CompressedMatrixSKS};

/// Element's stiffness matrix under static linear analysis.
pub trait StaticStiffness {
    fn k(&mut self) -> &CompressedMatrixSKS;
    fn k_printr(&self, n_exp: Dtype);
    fn k_string(&self, n_exp: Dtype) -> String;
}

/// Strain and stress data
pub trait SData {
    fn elem_size(&self) -> Dtype;
    fn nodes_ids(&self) -> Vec<usize>;
    fn strain_at_nodes(&mut self) -> Vec<Dtype>;
    fn stress_at_nodes(&mut self) -> Vec<Dtype>;
}

/// Export trait generate Part's informations and write it into files
/// Export trait generate txt or vtk files for output
pub trait Export {
    fn txt_writer(
        &self,
        target_file: &str,
        calc_time: std::time::Duration,
        n_exp: Dtype,
        energy: (Dtype, Dtype, Dtype),
    ) -> std::io::Result<bool>;

    fn vtk_writer(&mut self, target_file: &str, elem_type: &str) -> std::io::Result<bool>;
}
