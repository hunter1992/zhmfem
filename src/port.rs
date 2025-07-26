use crate::dtty::{basic::Dtype, matrix::CompressedMatrix};

/// K trait generate element's stiffness matrix under linear analysis.
/// Output stress/strain vector at some point in element.
pub trait K {
    fn k(&mut self) -> &CompressedMatrix;
    fn k_printer(&self, n_exp: Dtype);
    fn k_string(&self, n_exp: Dtype) -> String;

    fn strain_at_intpt(&mut self) -> Vec<Vec<Dtype>>;
    fn stress_at_intpt(&mut self) -> Vec<Vec<Dtype>>;

    //fn strain_at_nodes(&mut self) -> Vec<Vec<Dtype>>;
    //fn stress_at_nodes(&mut self) -> Vec<Vec<Dtype>>;

    fn info(&self, n_exp: Dtype) -> String;
    fn id(&self) -> usize;
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
