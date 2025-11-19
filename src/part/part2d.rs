use crate::calc::LinearEqs;
use crate::dtty::{
    basic::{ADtype, Dtype},
    matrix::CompressedMatrixSKS,
};
use crate::node::Node2D;
use crate::port::{Export, K};
use crate::tool::compress_matrix_sks;
use std::fmt::Write as _;
use std::io::{BufWriter, Write};
use std::sync::atomic::Ordering;
use std::sync::Arc;
use std::thread;
use std::time::Instant;
use vtkio::model::*;

/// Three generic const: N for N_NODE, F for N_FREEDOM, M for N_NODE in single element
/// In general, the value of F in a two-dimensional plane is 2
pub struct Part2D<'part2d, Elem: K, const N: usize, const F: usize, const M: usize>
where
    [[Dtype; N * F]; N * F]: Sized,
    [[Dtype; M * F]; M * F]: Sized,
{
    pub id: usize,
    pub nodes: &'part2d [Node2D],
    pub elems: &'part2d mut [Elem],
    pub cplds: &'part2d [Vec<usize>], //cplds: coupled nodes index
    pub k_matrix: Option<CompressedMatrixSKS>,
    pub assembly_time_consuming: Option<std::time::Duration>,
}

impl<'part2d, Elem: K, const N: usize, const F: usize, const M: usize>
    Part2D<'part2d, Elem, N, F, M>
where
    [[Dtype; N * F]; N * F]: Sized,
    [[Dtype; M * F]; M * F]: Sized,
{
    /// Construct a Part2D
    /// Args: nodes --- Vec ptr to all 2D nodes
    ///       elems --- Slice of all element stiffness matrix pointers
    ///       cplds --- Slice of the ID combinations of Node2D within a single cell
    pub fn new(
        id: usize,
        nodes: &'part2d [Node2D],
        elems: &'part2d mut [Elem],
        cplds: &'part2d [Vec<usize>],
    ) -> Self {
        Part2D {
            id,
            nodes,
            elems,
            cplds,
            k_matrix: None,
            assembly_time_consuming: None,
        }
    }

    /// Get nodes coordinates
    pub fn nodes_coordinate(&self) -> [Dtype; N * F] {
        let mut coords: [Dtype; N * F] = [0.0; N * F];
        for (idx, node) in self.nodes.iter().enumerate() {
            coords[idx * F] = node.coords[0];
            coords[idx * F + 1] = node.coords[1];
        }
        coords
    }

    /// Get nodes coordinates for vtk files
    pub fn vtk_nodes_coordinate(&self) -> Vec<Dtype> {
        let mut coords: Vec<Dtype> = vec![0.0; N * 3];
        let data = self.nodes_coordinate();
        for idx in 0..N {
            coords[idx * 3] = data[idx * 2];
            coords[idx * 3 + 1] = data[idx * 2 + 1];
            coords[idx * 3 + 2] = 0.0 as Dtype;
        }
        coords
    }

    /// Get displacement of all nodes
    pub fn nodes_displacement(&self) -> [Dtype; N * F] {
        let mut data: [Dtype; N * F] = [0.0; N * F];
        for (idx, node) in self.nodes.iter().enumerate() {
            data[idx * F] = node.displs.borrow()[0];
            data[idx * F + 1] = node.displs.borrow()[1];
        }
        data
    }

    /// Get displacement of all nodes for vtk files
    pub fn vtk_nodes_displacement(&self) -> Vec<Dtype> {
        let mut displs: Vec<Dtype> = vec![0.0; N * 3];
        let data = self.nodes_displacement();
        for idx in 0..N {
            displs[idx * 3] = data[idx * 2];
            displs[idx * 3 + 1] = data[idx * 2 + 1];
            displs[idx * 3 + 2] = 0.0 as Dtype;
        }
        displs
    }

    /// Get force of all nodes
    pub fn nodes_force(&self) -> [Dtype; N * F] {
        let mut data: [Dtype; N * F] = [0.0; N * F];
        for (idx, node) in self.nodes.iter().enumerate() {
            data[idx * F] = node.forces.borrow()[0];
            data[idx * F + 1] = node.forces.borrow()[1];
        }
        data
    }

    /// Get force of all nodes for vtk files
    pub fn vtk_nodes_force(&self) -> Vec<Dtype> {
        let mut forces: Vec<Dtype> = vec![0.0; N * 3];
        let data = self.nodes_force();
        for idx in 0..N {
            forces[idx * 3] = data[idx * 2];
            forces[idx * 3 + 1] = data[idx * 2 + 1];
            forces[idx * 3 + 2] = 0.0 as Dtype;
        }
        forces
    }

    /// Calculate part's global stiffness matrix
    /// The value of arg parallel_or_singllel can be:
    /// "parallel" or "p" or "singllel" or "s"
    pub fn k(&mut self, parallel_or_singllel: &str, cpu_cores: usize) -> &CompressedMatrixSKS {
        if self.k_matrix.is_none() {
            if self.cplds.len() != self.elems.len() {
                println!("\n---> Error! From Part2D.k func.");
                println!("     The count of elements not eq to K mat size.");
                panic!("---> Assembly global K failed!");
            }

            match parallel_or_singllel {
                "s" => {
                    self.k_matrix = Some(self.k_singllel(cpu_cores));
                }
                "p" => self.k_matrix = Some(self.k_parallel(cpu_cores)),
                "singllel" => self.k_matrix = Some(self.k_singllel(cpu_cores)),
                "parallel" => self.k_matrix = Some(self.k_parallel(cpu_cores)),
                _ => {
                    if self.cplds.len() >= 1000 {
                        self.k_matrix = Some(self.k_parallel(cpu_cores))
                    } else {
                        self.k_matrix = Some(self.k_singllel(cpu_cores))
                    }
                }
            }
        }
        self.k_matrix.as_ref().unwrap()
    }

    fn k_singllel(&mut self, _cpu_cores: usize) -> CompressedMatrixSKS {
        println!(
            "\n>>> Assembling Part2D(#{})'s stiffness matrix K{}(size: {}x{}) in single thread ......",
            self.id, self.id, N*F, N*F 
        );

        let elems: Vec<_> = self.elems.iter_mut().map(|elem| elem.k()).collect();

        let mut part_stiffness_mat: [[Dtype; N * F]; N * F] = [[0.0; N * F]; N * F];

        let timing_start = Instant::now();

        self.cplds
            .iter()
            .zip(elems.iter())
            .map(|(cpld, &elem)| {
                let elem_stiffness_mat = elem.recover::<{ M * F }>();
                for idx in 0..M {
                    for idy in 0..M {
                        for row in 0..F {
                            for col in 0..F {
                                part_stiffness_mat[cpld[idx] * F + row][cpld[idy] * F + col] +=
                                    elem_stiffness_mat[idx * F + row][idy * F + col];
                            }
                        }
                    }
                }
            })
            .count();

        let time_consuming = timing_start.elapsed();
        println!(
            "\n>>> Assembly Down!\n        time consuming: {:?}",
            time_consuming
        );
        self.assembly_time_consuming = Some(time_consuming);

        compress_matrix_sks(&part_stiffness_mat)
    }

    fn k_parallel(&mut self, cpu_cores: usize) -> CompressedMatrixSKS {
        println!(
            "\n>>> Assembling Part2D(#{})'s stiffness matrix K{} (size: {}x{}) in {} threads ......",
            self.id, self.id, N*F, N*F, cpu_cores
        );

        let elems: Vec<_> = self.elems.iter_mut().map(|elem| elem.k()).collect();
        let mut part_stiffness_mat: [[Dtype; N * F]; N * F] = [[0.0; N * F]; N * F];
        let mat: Arc<Vec<Vec<ADtype>>> = Arc::new(
            (0..F * N)
                .map(|_| (0..F * N).map(|_| ADtype::new(0.0)).collect())
                .collect(),
        );

        let density: usize = (((self.cplds.len() as f64) / (cpu_cores as f64)).ceil()) as usize;
        let elemss: Vec<_> = elems.chunks(density).collect();
        let cpldss: Vec<_> = self.cplds.chunks(density).collect();

        let timing_start = Instant::now();

        for n in 0..cpu_cores {
            let elems = elemss[n];
            let cplds = cpldss[n];
            let mat_clone = Arc::clone(&mat);
            thread::scope(|s| {
                s.spawn(|| {
                    cplds
                        .iter()
                        .zip(elems.iter())
                        .map(|(cpld, &elem)| {
                            for idx in 0..M {
                                for idy in 0..M {
                                    for row in 0..F {
                                        for col in 0..F {
                                            let elem_stiffness_mat = elem.recover::<{ M * F }>();
                                            mat_clone[cpld[idx] * F + row][cpld[idy] * F + col]
                                                .fetch_add(
                                                    elem_stiffness_mat[idx * F + row]
                                                        [idy * F + col],
                                                    Ordering::Relaxed,
                                                );
                                        }
                                    }
                                }
                            }
                        })
                        .count();
                });
            });
        }

        let time_consuming = timing_start.elapsed();
        println!(
            "\n>>> Assembly Down!\n        time consuming: {:?}",
            time_consuming
        );
        self.assembly_time_consuming = Some(time_consuming);
        for x in 0..F * N {
            for y in 0..F * N {
                part_stiffness_mat[x][y] = mat[x][y].load(Ordering::Relaxed);
            }
        }
        compress_matrix_sks(&part_stiffness_mat)
    }

    /// Print part's global stiffness matrix
    pub fn k_printer(&mut self, parallel_or_singllel: &str, cpu_cores: usize, n_exp: Dtype) {
        if self.k_matrix.is_none() {
            self.k(parallel_or_singllel, cpu_cores);
        }

        let k_matrix = self.k_matrix.clone().unwrap().recover::<{ N * F }>();
        print!("\nPart #{}  K =  (* 10^{})\n[", self.id, n_exp as u8);
        for row in 0..(N * F) {
            if row == 0 {
                print!("[");
            } else {
                print!(" [");
            }
            for col in 0..(N * F) {
                print!(
                    " {:>-13.6} ",
                    k_matrix[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
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

    /// Return Part2D's stiffness matrix's formeted string
    pub fn k_string(&self, n_exp: Dtype) -> String {
        let mut k_matrix = String::new();
        write!(
            k_matrix,
            "\nPart #{}  K =  (* 10^{})\n",
            self.id, n_exp as u8
        )
        .expect("Write Part2D Stiffness Matrix Failed!");

        let k: [[Dtype; N * F]; N * F] = *self.k_matrix.clone().unwrap().recover();
        for row in 0..(N * F) {
            if row == 0 {
                write!(k_matrix, "[[").expect("Write Part2D Stiffness Matrix Failed!");
            } else {
                write!(k_matrix, " [").expect("Write Part2D Stiffness Matrix Failed!");
            }
            for col in 0..(N * F) {
                write!(
                    k_matrix,
                    " {:>-13.6} ",
                    k[row][col] / (10.0_f64.powf(n_exp as f64)) as Dtype
                )
                .expect("Write Part2D Stiffness Matrix Failed!");
            }
            if row == { N * F - 1 } {
                write!(k_matrix, "]]").expect("Write Part2D Stiffness Matrix Failed!");
            } else {
                write!(k_matrix, " ]\n").expect("Write Part2D Stiffness Matrix Failed!");
            }
        }
        write!(k_matrix, "\n").expect("Write Part2D Stiffness Matrix Failed!");
        k_matrix
    }

    /// Write the disp and force result into nodes
    pub fn write_result(&self, slv: &LinearEqs<{ N * F }>) {
        let disp = slv.disps;
        let force = slv.external_force.unwrap();
        for (idx, node) in self.nodes.iter().enumerate() {
            node.displs.borrow_mut()[0] = disp[idx * 2];
            node.displs.borrow_mut()[1] = disp[idx * 2 + 1];
            node.forces.borrow_mut()[0] = force[idx * 2];
            node.forces.borrow_mut()[1] = force[idx * 2 + 1];
        }
    }

    /// Returns the strain value for each constant strain element
    pub fn elem_strain(&mut self) -> Vec<Vec<Dtype>> {
        let mut strain: Vec<Vec<Dtype>> = Vec::with_capacity(self.cplds.len());
        self.elems
            .iter_mut()
            .map(|elem| strain.push(elem.strain_at_intpt().remove(0)))
            .count();
        strain
    }

    /// Returns the stress value for each constant strain element
    pub fn elem_stress(&mut self) -> Vec<Vec<Dtype>> {
        let mut stress: Vec<Vec<Dtype>> = Vec::with_capacity(self.cplds.len());
        self.elems
            .iter_mut()
            .map(|elem| stress.push(elem.stress_at_intpt().remove(0)))
            .count();
        stress
    }

    /// Wtrie triangle elements into vtk file
    pub fn vtk_tri2d3n_legacy(&mut self, vtk_file_name: &str) -> Vtk {
        let elem_num: usize = self.cplds.len();
        let coords = self.vtk_nodes_coordinate();
        let displs = self.vtk_nodes_displacement();
        let forces = self.vtk_nodes_force();

        let mut vtk_cpld: Vec<u32> = vec![0; 4 * elem_num];
        for idx in 0..elem_num {
            vtk_cpld[idx * 4] = 3;
            vtk_cpld[idx * 4 + 1] = self.cplds[idx][0] as u32;
            vtk_cpld[idx * 4 + 2] = self.cplds[idx][1] as u32;
            vtk_cpld[idx * 4 + 3] = self.cplds[idx][2] as u32;
        }

        let strain = self.elem_strain();
        let stress = self.elem_stress();

        let mut cell_e11: Vec<Dtype> = vec![0.0; elem_num];
        let mut cell_e22: Vec<Dtype> = vec![0.0; elem_num];
        let mut cell_e12: Vec<Dtype> = vec![0.0; elem_num];

        let mut cell_s11: Vec<Dtype> = vec![0.0; elem_num];
        let mut cell_s22: Vec<Dtype> = vec![0.0; elem_num];
        let mut cell_s12: Vec<Dtype> = vec![0.0; elem_num];

        for idx in 0..elem_num {
            cell_e11[idx] = strain[idx][0];
            cell_e22[idx] = strain[idx][1];
            cell_e12[idx] = strain[idx][2];
            cell_s11[idx] = stress[idx][0];
            cell_s22[idx] = stress[idx][1];
            cell_s12[idx] = stress[idx][2];
        }

        Vtk {
            version: Version { major: 4, minor: 2 },
            title: String::from(vtk_file_name),
            byte_order: ByteOrder::BigEndian,
            file_path: None,
            data: DataSet::inline(UnstructuredGridPiece {
                // points: IOBuffer::F32(coords),
                points: IOBuffer::F64(coords),
                cells: Cells {
                    cell_verts: VertexNumbers::Legacy {
                        num_cells: elem_num as u32,
                        vertices: vtk_cpld,
                    },
                    types: vec![CellType::Triangle; elem_num],
                },
                data: Attributes {
                    point: vec![
                        Attribute::vectors("displacement").with_data(displs),
                        Attribute::vectors("force").with_data(forces),
                    ],
                    cell: vec![
                        Attribute::scalars("e_xx", 1).with_data(cell_e11),
                        Attribute::scalars("e_yy", 1).with_data(cell_e22),
                        Attribute::scalars("e_xy", 1).with_data(cell_e12),
                        Attribute::scalars("s_xx", 1).with_data(cell_s11),
                        Attribute::scalars("s_yy", 1).with_data(cell_s22),
                        Attribute::scalars("s_xy", 1).with_data(cell_s12),
                    ],
                },
            }),
        }
    }

    pub fn vtk_quad2d4n_legacy(&mut self, vtk_file_name: &str) -> Vtk {
        let elem_num: usize = self.cplds.len();
        let coords = self.vtk_nodes_coordinate();
        let displs = self.vtk_nodes_displacement();
        let forces = self.vtk_nodes_force();

        let mut vtk_cpld: Vec<u32> = vec![0; 5 * elem_num];
        for idx in 0..elem_num {
            vtk_cpld[idx * 5] = 4;
            vtk_cpld[idx * 5 + 1] = self.cplds[idx][0] as u32;
            vtk_cpld[idx * 5 + 2] = self.cplds[idx][1] as u32;
            vtk_cpld[idx * 5 + 3] = self.cplds[idx][2] as u32;
            vtk_cpld[idx * 5 + 4] = self.cplds[idx][3] as u32;
        }

        Vtk {
            version: Version { major: 4, minor: 2 },
            title: String::from(vtk_file_name),
            byte_order: ByteOrder::BigEndian,
            file_path: None,
            data: DataSet::inline(UnstructuredGridPiece {
                // points: IOBuffer::F32(coords),
                points: IOBuffer::F64(coords),
                cells: Cells {
                    cell_verts: VertexNumbers::Legacy {
                        num_cells: elem_num as u32,
                        vertices: vtk_cpld,
                    },
                    types: vec![CellType::Quad; elem_num],
                },
                data: Attributes {
                    point: vec![
                        Attribute::vectors("displacement").with_data(displs),
                        Attribute::vectors("force").with_data(forces),
                    ],
                    cell: vec![],
                },
            }),
        }
    }
}
impl<'part2d, Elem: K, const N: usize, const F: usize, const M: usize> Export
    for Part2D<'part2d, Elem, N, F, M>
where
    [[Dtype; N * F]; N * F]: Sized,
    [[Dtype; M * F]; M * F]: Sized,
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
        write!(text_writer, "\n").expect("Write parts' result into txt file failed!!!");
        write!(text_writer, "\n>>> System energy:\n")
            .expect("Write parts' result into txt file failed!!!");
        write!(text_writer, "\t\t\tE_d: {:-13.6} (deform energy)\n", deform)
            .expect("Write parts' result into txt file failed!!!");
        write!(text_writer, "\t\t\tW_f: {:-13.6} (exforce works)\n", exwork)
            .expect("Write parts' result into txt file failed!!!");
        write!(
            text_writer,
            "\t\t\tE_p: {:-13.6} (potential energy)\n",
            potential
        )
        .expect("Write parts' result into txt file failed!!!");

        write!(
            text_writer,
            "\n>>> Solver time consuming: {:?}\n",
            calc_time
        )
        .expect("Write txt file error!");

        write!(text_writer, "\n>>> Part2D stiffness matrix:\n").expect("Write info failed!");
        write!(text_writer, "{}", self.k_string(n_exp)).expect("Write info failed!");

        write!(text_writer, "\n>>> Details of each element:")
            .expect("Write parts' result into txt file failed!!!");
        for elem in self.elems.iter() {
            write!(text_writer, "{}\n", elem.info(n_exp)).expect("Write info failed!");
        }
        text_writer.flush().expect("!!! Flush txt file failed!");
        println!("    Down!");
        Ok(true)
    }

    fn vtk_writer(&mut self, target_file: &str, elem_type: &str) -> std::io::Result<bool> {
        let vtk_file = std::fs::File::create(target_file).unwrap();
        let mut vtk_writer = BufWriter::new(vtk_file);
        let mut vtk_string = String::new();

        match elem_type {
            "Tri2D3N" => {
                println!("\n>>> Writing calc results into vtk file ......");
                let vtk_triangle = self.vtk_tri2d3n_legacy(elem_type);
                vtk_triangle.write_legacy_ascii(&mut vtk_string).expect(
                    "!!! Construct VTK string for Tri2D3N failed! from zhmfem/src/part/part2d.rs",
                );

                write!(vtk_writer, "{}", vtk_string.as_str())
                    .expect("!!! Write .vtk file failed! from zhmfem/src/part/part2d.rs");
                println!("    Down!");

                Ok(true)
            }

            "Quad2D4N" => {
                println!("\n>>> Writing calc results into vtk file ......");
                let vtk_quadrilateral = self.vtk_quad2d4n_legacy(elem_type);
                vtk_quadrilateral
                    .write_legacy_ascii(&mut vtk_string)
                    .expect(
                    "!!! Construct VTK string for Tri2D3N failed! from zhmfem/src/part/part2d.rs",
                );

                write!(vtk_writer, "{}", vtk_string.as_str())
                    .expect("!!! Write .vtk file failed! from zhmfem/src/part/part2d.rs");
                println!("    Down!");

                Ok(true)
            }
            _ => Ok(false),
        }
    }
}
